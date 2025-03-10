/*
 * mcmc_parallel.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include <atomic>
#include "mcmc_build.h"

/*
struct LoglikeCalculator : public RcppParallel::Worker
{
    MCMC &mcmc;
    std::atomic<bool> &has_error;
    std::string &error_message;
    double new_ll;

    LoglikeCalculator(MCMC &mcmc, std::atomic<bool> &has_error, std::string &error_message)
        : mcmc(mcmc), has_error(has_error), error_message(error_message), new_ll(0.0) {}

    LoglikeCalculator(const LoglikeCalculator &other, RcppParallel::Split)
        : mcmc(other.mcmc), has_error(other.has_error), error_message(other.error_message), new_ll(0.0) {}

    void operator()(size_t begin, size_t end)
    {
        for (size_t H = begin; H < end; ++H)
        {
            if (has_error.load(std::memory_order_relaxed))
                return;

            try
            {
                Household &house = mcmc.get_household(H);

                // Update house loglikelihood
                house.set_loglike(mcmc.compute_house_loglike(house));
                new_ll += house.get_loglike();

                if (mcmc.is_augm_model())
                    new_ll += mcmc.compute_augm_model_loglike(house);
            }
            catch (const std::exception &e)
            {
                if (!has_error.exchange(true))
                { // Only set error once
                    error_message = "While computing likelihood for household " +
                                    mcmc.get_household(H).house_id() + ": " + e.what();
                }
            }
        }
    }

    void join(const LoglikeCalculator &other)
    {
        new_ll += other.new_ll;
    }
};

void MCMC::update_simple_global_loglike()
{
    std::atomic<bool> has_error(false);
    std::string error_message;
    double new_ll = 0.0;

    LoglikeCalculator calc(*this, has_error, error_message);
    RcppParallel::parallelReduce(0, n_houses(), calc);

    if (has_error.load())
        tools::error(error_message);

    pImpl->set_global_loglike(calc.new_ll);
}
*/

/* Set the number of threads for parrell computations (default is 1). */
void MCMC::Impl::set_n_threads()
{
    if (setup_n_threads == 0)
        throw std::invalid_argument("N. of threads must be greater than 0.");

#ifdef _OPENMP
    unsigned int max_threads = omp_get_max_threads();
    if (setup_n_threads > max_threads)
    {
        tools::warning(
            "N. of threads reduced to " + std::to_string(max_threads) + " available threads.");
        setup_n_threads = max_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(setup_n_threads);

#else
    if (setup_n_threads != 1)
        tools::warning("Parallelization with OpenMP not supported.");

    setup_n_threads = 1;
#endif
}

/* Assign seed for random generators. */
void MCMC::Impl::set_prng_seed()
{
    // seed for the unique thread rand num generator
    prng.seed(setup_seed);
    _func_prng_ = &Impl::_prng_unique_thread_;

#ifdef _OPENMP
    // seeds for the case of parallelization
    if (setup_n_threads > 1)
    {
        // create different generators for different parallel threads
        // note: they are used only to parallelize data augmentantion
        _func_prng_ = &Impl::_prng_parallel_threads_;

        std::uniform_int_distribution<unsigned long> seed_generator(0, 0XFFFFFFFF);

        prng_parallel.reserve(setup_n_threads);
        for (unsigned n = 0; n < setup_n_threads; ++n)
            prng_parallel[n].seed(seed_generator(prng));
    }
#endif
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Compute new loglikelihood, store old value in buffer. Same for households. */
void MCMC::update_simple_global_loglike()
{
    double new_ll = 0.0;
    std::atomic<bool> has_error(false);

#pragma omp parallel for reduction(+ : new_ll)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        try
        {
            Household &house = get_household(H);

            // update house loglikelihood
            house.set_loglike(compute_house_loglike(house));
            new_ll += house.get_loglike();
            if (is_augm_model())
                new_ll += compute_augm_model_loglike(house);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    _error_message_ = "While computing likelihood for household " + get_household(H).house_id() + ": " + e.what();
            }
        }
    }
    if (!_error_message_.empty())
        tools::error(_error_message_);

    pImpl->set_global_loglike(new_ll);
}

/* Store all buffers (also for extra dependecies) and update values of loglikelihood. */
void MCMC::update_extra_global_loglike()
{
    double new_ll = 0.0;
    std::atomic<bool> has_error(false);

#pragma omp parallel for reduction(+ : new_ll)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        try
        {
            Household &house = get_household(H);

            // buffer all extra dependencies
            house.buffer_all_extra();
            for (const indiv &i : house.infected_status_list())
                update_extra_dependecies(house, i);

            // update house loglikelihood
            house.set_loglike(compute_house_loglike(house));
            new_ll += house.get_loglike();
            if (is_augm_model())
                new_ll += compute_augm_model_loglike(house);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    _error_message_ = "While computing likelihood for household " + get_household(H).house_id() + " : " + e.what();
            }
        }
    }
    if (!_error_message_.empty())
        tools::error(_error_message_);

    pImpl->set_global_loglike(new_ll);
}

/* Move augmented variables from all households. */
void MCMC::move_all_augm()
{
    double delta_ll = 0.;
    std::atomic<bool> has_error(false);

#pragma omp parallel for reduction(+ : delta_ll)
    for (size_t H = 0; H < n_houses(); ++H)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        try
        {
            delta_ll += move_augm_single_house(get_household(H));
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    _error_message_ = "While moving augmented variables in household " + get_household(H).house_id() + ": " + e.what();
            }
        }
    }

    if (!_error_message_.empty())
        tools::error(_error_message_);

    pImpl->increment_global_loglike(delta_ll);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Perform simulation. */
void MCMC::aux_simulate()
{
    // FIXME : memory usage is suboptimal:
    // at a given moment 3 memory copies of dataset-like objects
    std::vector<Household> hhvec = {};

    try
    {
        std::atomic<bool> has_error(false);

        //==========================//
        //===== Initialization =====//
        //==========================//

        pImpl->assert_proper_dataset(); // check if dataset is empty
        pImpl->setup_n_chains = 0;      // use 0 chains as a state variable

        pImpl->setup_n_threads = 1; // impose no parallelization
        pImpl->set_n_threads();

        pImpl->setup_rand_init_params = false; // disable random parameter initialization
        pImpl->init_parameters();              // check no undefined parameters
        init_infec_profile();                  // initialize user infectivity profile
        init_augmented_variables();            // initialize augmented dependencies
        init_model_for_simulation();           // initialize model for simulation

        //======================//
        //===== Simulation =====//
        //======================//

        progress::Bar progressbar(n_houses(), "simulate", pImpl->setup_progress);

        for (size_t H = 0; H < n_houses(); ++H)
        {
            if (has_error.load(std::memory_order_relaxed))
                continue;

            try
            {
                Household house_copy = get_household(H);
                simulate_single_house(house_copy);
                
#ifdef _HOUSEPI_RCPP
                pImpl->_aux_export_dataset_update_(house_copy);
#else
                pImpl->_aux_print_dataset_update_(house_copy);
#endif
            }
            catch (const std::exception &e)
            {
                if (!has_error.exchange(true))
                    _error_message_ = "During simulation for household " + get_household(H).house_id() + " : " + e.what();
            }
            progressbar.update(H);
        }

        if (!_error_message_.empty())
            throw std::runtime_error(_error_message_);

        progressbar.finalize();
    }
    catch (std::exception &e)
    {
        tools::error(e.what());
    }
}
