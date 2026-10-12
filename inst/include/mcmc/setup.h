/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <atomic>
#include <cstddef>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/random.h"

//////////////////////////////  Setup  //////////////////////////////

namespace housepi
{
    /// Split workers between the two levels of parallel dispatch.
    struct Budget
    {
        std::size_t outer = 1; /// Workers spread across chains.
        std::size_t inner = 1; /// Workers spread within one chain.
    };

    /// Configure an MCMC run.
    ///
    /// * holds user, advanced and derived options;
    /// * owns the random number generators and the iteration count of a run.
    struct Setup
    {
        //////  User config  //////

        std::int32_t seed;                    /// Seed of the pseudo-random number generators.
        unsigned n_threads;                   /// Number of threads for parallel computations.
        unsigned long n_steps;                /// Number of iterations per chain after burn-in, a multiple of n_thin.
        unsigned long n_burn;                 /// Number of burn-in iterations to discard, run after adaptation.
        unsigned long n_thin;                 /// Thinning interval, storing every n_thin-th iteration.
        unsigned long n_adapt;                /// Number of adaptive MCMC steps, 0 for no adaptation.
        unsigned long n_chains;               /// Number of independent MCMC chains.
        std::vector<std::string> id_draws;    /// Individual IDs whose latent draws are traced.
        std::vector<std::string> house_draws; /// Household IDs whose infected members are traced.

        //////  Advanced config  //////

        unsigned long n_adapt_window; /// Effective window for the smoothed accept-rate.
        double time_step;             /// Real-time integration step for simulation.

        //////  Derived config  //////

        Budget budget;               /// Workers allotted to the chain and the per-chain dispatch level.
        unsigned long c_sample_size; /// Sample size for a single chain after thinning.
        double adapt_smoothing;      /// Smoothing factor of the accept-rate average, 1 / n_adapt_window.

        //////  Run config  //////

        const std::atomic<bool> *stop_flag = nullptr; /// Flag shared by the engines of one run, null outside it.
        std::atomic<unsigned long> n_iter{0};         /// Iterations counted in the run, written by one thread.

        //////  Constructor  //////

        Setup();

        //////  Setter  //////

        /// Set the number of threads, capped at those available.
        void set_n_threads();
        /// Reset every field to its default.
        void reset_defaults();
        /// Build setup from input table.
        void build(const housepi::str_table &input);
        /// Derive the sample size, then settle the adaptation steps.
        /// * nothing adaptive -> adaptation steps turn into burn-in when `burn` is 0, else are dropped.
        /// @note Counts are asserted R-side, `steps` extended to a multiple of `thin`.
        void init_mcmc_configuration(bool is_adaptive);

        //////  Random  //////

        /// Seed one generator per thread.
        ///
        /// * thread 0 takes `seed`;
        /// * each next thread jumps 2^128 steps ahead of the previous one.
        void set_prng_seed();
        /// Return generator of the calling thread.
        housepi::rnd_gen &prng() const;

        /// Sample a uniform-distributed variable in [0, 1).
        double uniform() const { return housepi::standard_uniform(prng()); }
        /// Sample an exponential-distributed variable with rate 1.
        double exponential() const { return housepi::standard_exponential(prng()); }

        //////  Progress  //////

        /// Count one iteration.
        /// * stop flag raised -> throws parallel::Cancelled.
        void count_iter();
        /// Return fraction of the run's iterations counted, from 0 to 1.
        double fraction_iter() const;

    private:
        mutable std::vector<housepi::rnd_gen> prng_parallel_ = {}; /// One generator per thread.

        /// Associate non-default user setup.
        void helper_non_default_build_(const housepi::str_table &input);
    };
}
