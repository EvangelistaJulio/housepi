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

#include "tools/error.h"
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : Impl : simulation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _HOUSEPI_RCPP
/// Perform a unique repetition.
Rcpp::DataFrame MCMC::Impl::simulate_repetition()
{
    using namespace Rcpp;

    const size_t n_rows = dataset.n_individuals();

    IntegerVector df_status = IntegerVector(n_rows);

    std::vector<NumericVector> df_augm;
    df_augm.reserve(augm.list_export.size());
    for (size_t c = 0; c < augm.list_export.size(); ++c)
        df_augm.emplace_back(n_rows);

    size_t r = 0;
    for (Household house : dataset) // NOTE: deep copy of household
        try
        {
            simulate_single_household_(house); // inplace simulation
            for (const Individual &member : house)
            {
                df_status[r] = member.status();
                for (size_t c = 0; c < augm.list_export.size(); ++c)
                    df_augm[c][r] = hsp_to_rcpp(member.augm(augm.list_export[c]));
                ++r;
            }
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While simulating " + house.err_name() + ":\n" +
                std::string(e.what()));
        }

    DataFrame output = DataFrame::create(_["status"] = df_status);
    for (size_t c = 0; c < augm.list_export.size(); ++c)
        output.push_back(df_augm[c], augm.list_export[c]);

    return output;
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
        hsp::error(e.what());
    }

    return List::create();
}

#endif
