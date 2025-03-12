/*
 * mcmc_parallel.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include <atomic>
#include "mcmc_build.h"

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
    for (ulong H = 0; H < n_houses(); ++H)
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
    for (ulong H = 0; H < n_houses(); ++H)
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
    for (ulong H = 0; H < n_houses(); ++H)
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
