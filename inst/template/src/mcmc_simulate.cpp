/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_build_param.cpp
 *
 */

#include <atomic>

#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : Impl : simulation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Auxiliary method to initialize simulations.
// UPDATE: once I allow for multi-treading and a flag than n_chains = 0, I can drop this...
void MCMC::Impl::aux_init_simulate_()
{
    assert_proper_dataset_();       // check if dataset is empty
    setup.n_chains = 0;             // use 0 chains as a state variable
    setup.n_threads = 1;            // impose no parallelization
    setup.set_n_threads();          // change accordingly
    init_model_for_simulation_();   // initialize model for simulation
}

#ifdef _HOUSEPI_RCPP
/// Perform a unique repetition.
Rcpp::DataFrame MCMC::Impl::simulate_repetition()
{
    using namespace Rcpp;

    const size_t n_cols = dataset.n_individuals();

    int_vec df_status = int_vec(n_cols);
    std::vector<dbl_vec> df_augm(augm.list_export.size(), dbl_vec(n_cols));

    size_t row = 0;
    for (Household house : dataset) // copies of households
        try
        {
            simulate_single_household_(house); // inplace simulation
            for (const Individual &member : house)
            {
                df_status[row] = member.status();
                for (size_t i = 0; i < augm.list_export.size(); ++i)
                    df_augm[i][row] = member.augm(augm.list_export[i]);
                ++row;
            }
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While simulating " + house.err_name() + ":\n" + 
                std::string(e.what()));
        }

    List output = List::create(_["status"] = wrap(std::move(df_status)));
    for (size_t i = 0; i < augm.list_export.size(); ++i)
        output[augm.list_export[i]] = wrap(std::move(df_augm[i]));
    return DataFrame(output);
}

/// Perform simulation and return simulated dataset.
Rcpp::List MCMC::Impl::simulate(size_t n_rep)
{
    using namespace Rcpp;
    std::atomic<bool> has_error(false);

    try
    {
        // initialization
        aux_init_simulate_();
        setup.prog_open(dataset.n_houses() * n_rep);

        // run simulations
        std::vector<DataFrame> results(n_rep);

#pragma omp parallel for
        for (size_t rep = 0; rep < n_rep; ++rep)
        {
            if (has_error.load(std::memory_order_relaxed))
                continue;
            try
            {
                results[rep] = simulate_repetition();
            }
            catch (const std::exception &e)
            {
#pragma omp critical
                {
                    if (!has_error.exchange(true))
                        error_message = e.what();
                }
            }
        }
        if (!error_message.empty())
            throw std::runtime_error(error_message);

        // convert to Rcpp List
        List output(n_rep);
        for (size_t rep = 0; rep < n_rep; ++rep)
            output[rep] = std::move(results[rep]);

        // close progress bar
        setup.prog_close();

        return output;
    }
    catch (std::exception &e)
    {
        tools::error(e.what());
    }

    return List::create();
}

#else
/// Print the simulated dataset to the provided file.
void MCMC::Impl::print_simulated_data_wide(const std::string &pathtofile)
{
    try
    {
        // initialization
        setup.prog_declare("simulate");
        dataset.print_open(pathtofile);
        aux_init_simulate_();

        setup.prog_open(dataset.n_houses()); // open progress bar

        // run simulation
        // NOTE: deep copy of the household
        for (Household house : dataset)
            try
            {
                simulate_single_household_(house); // inplace simulation
                dataset.print_update(house);       // print the simulated dataset
                setup.prog_update();               // update progress bar
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While simulating " + house.err_name() + ":\n" +
                    e.what());
            }

        setup.prog_close();    // close progress bar
        dataset.print_close(); // finalize output
    }
    catch (const std::exception &e)
    {
        setup.prog_crash();
        tools::error(e.what());
    }
}
#endif
