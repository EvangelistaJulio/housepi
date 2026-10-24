/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/parallel.h"
#include "tools/random.h"

//////////////////////////////  Settings  //////////////////////////////

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
    /// * owns the random number generators and the step count of a run.
    struct Settings
    {
        //////  Types  //////

        using Step = std::uint64_t;  /// Count or index of MCMC steps.
        using Chain = std::uint64_t; /// Count or index of MCMC chains.

        //////  User config  //////

        std::int32_t seed;                    /// Seed of the pseudo-random number generators.
        unsigned n_threads;                   /// Number of threads for parallel computations.
        Step n_steps;                         /// Number of MCMC steps per chain after the warm-up, a multiple of n_thin.
        Step n_thin;                          /// Thinning interval, storing every n_thin-th step.
        Step n_adapt;                         /// Number of adaptive MCMC steps, 0 for no adaptation.
        Chain n_chains;                       /// Number of independent MCMC chains.
        std::vector<std::string> id_draws;    /// Individual IDs whose latent draws are traced.
        std::vector<std::string> house_draws; /// Household IDs whose infected members are traced.
        double time_step;                     /// Real-time integration step for simulation.

        //////  Derived config  //////

        Budget budget;      /// Workers allotted to the chain and the per-chain dispatch level.
        Step c_sample_size; /// Sample size for a single chain after thinning.
        Step n_burn;        /// Number of burn-in steps, the adaptation budget when nothing adapts.

        //////  Run config  //////

        const std::atomic<bool> *stop_flag = nullptr; /// Flag shared by the engines of one run, null outside it.
        std::atomic<Step> n_done{0};                  /// Steps counted in the run, written by one thread.

        //////  Constructor  //////

        Settings() { reset_defaults(); }

        //////  Setter  //////

        /// Set the number of threads, capped at those available.
        void set_n_threads();
        /// Reset every field to its default.
        void reset_defaults();
        /// Build settings from input table.
        void build(const housepi::str_table &input);
        /// Derive the sample size, then settle the adaptation steps.
        /// * nothing adaptive -> adaptation steps turn into burn-in.
        /// @note Counts are asserted R-side, `steps` extended to a multiple of `thin`.
        void init_mcmc_configuration(bool is_adaptive);

        //////  Random  //////

        /// Seed one generator per thread.
        ///
        /// * thread 0 takes `seed`;
        /// * each next thread jumps 2^128 steps ahead of the previous one.
        void set_prng_seed();
        /// Return slot of the calling thread in per-thread pools, below n_threads.
        /// * one worker within a chain -> 0;
        /// * otherwise -> worker index in the arena.
        std::size_t thread_slot() const { return budget.inner == 1 ? 0 : housepi::parallel::current_thread_index(); }
        /// Return generator of the calling thread.
        housepi::RndGen &prng() const { return prng_parallel_[thread_slot()]; }

        /// Sample a uniform-distributed variable in [0, 1).
        double uniform() const { return housepi::standard_uniform(prng()); }
        /// Sample an exponential-distributed variable with rate 1.
        double exponential() const { return housepi::standard_exponential(prng()); }

        //////  Progress  //////

        /// Count one step.
        /// * stop flag raised -> throws parallel::Cancelled.
        void count_step()
        {
            if (stop_flag && stop_flag->load(std::memory_order_relaxed))
                throw parallel::Cancelled{};
            n_done.store(n_done.load(std::memory_order_relaxed) + 1, std::memory_order_relaxed);
        }
        /// Return fraction of the run's steps counted, from 0 to 1.
        double fraction_done() const;

    private:
        mutable std::vector<housepi::RndGen> prng_parallel_ = {}; /// One generator per thread.

        /// Associate non-default user settings.
        void helper_non_default_build_(const housepi::str_table &input);
    };
}
