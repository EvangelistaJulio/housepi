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
        /// Workers spread across chains.
        std::size_t outer = 1;
        /// Workers spread within one chain.
        std::size_t inner = 1;
    };

    /// Configure an MCMC run.
    ///
    /// * holds user, advanced and derived options;
    /// * owns the random number generators and the step count of a run.
    struct Settings
    {
        //////  Types  //////

        /// Count or index of MCMC steps.
        using Step = std::uint64_t;
        /// Count or index of MCMC chains.
        using Chain = std::uint64_t;

        //////  User config  //////

        std::int32_t seed;
        unsigned n_threads;
        /// Number of MCMC steps per chain after the warm-up, a multiple of n_thin.
        Step n_steps;
        /// Thinning interval, storing every n_thin-th step.
        Step n_thin;
        /// Number of adaptive MCMC steps, 0 for no adaptation.
        Step n_adapt;
        Chain n_chains;
        /// Individual IDs whose latent draws are traced.
        std::vector<std::string> id_draws;
        /// Household IDs whose infected members are traced.
        std::vector<std::string> house_draws;
        /// Real-time integration step for simulation.
        double time_step;

        //////  Derived config  //////

        /// Workers allotted to the chain and the per-chain dispatch level.
        Budget budget;
        /// Sample size for a single chain after thinning.
        Step c_sample_size;
        /// Number of burn-in steps, the adaptation budget when nothing adapts.
        Step n_burn;

        //////  Run config  //////

        /// Flag shared by the engines and chains of one run, null outside it.
        std::atomic<bool> *stop_flag = nullptr;
        /// Steps counted in the run, by every chain.
        std::atomic<Step> n_done{0};

        //////  Constructor  //////

        Settings() { reset_defaults(); }

        //////  Setter  //////

        /// Set the number of threads, capped at those available.
        void set_n_threads();
        /// Reset every field to its default.
        void reset_defaults();
        /// Build settings from input table.
        void build(const housepi::str_table &input);
        /// Split threads, then derive the sample size and adaptation steps.
        /// * threads -> `gcd(chains, threads)` chains at once, the rest to the households of each;
        /// * nothing adaptive -> adaptation steps turn into burn-in.
        /// @note Counts are asserted R-side, `steps` extended to a multiple of `thin`.
        void init_mcmc_configuration(bool is_adaptive);

        //////  Random  //////

        /// Seed one stream per chain, allot one household generator per thread.
        void set_prng_seed();
        /// Return slot of the calling thread in per-thread pools, below `budget.inner`.
        /// * one worker within a chain -> 0;
        /// * otherwise -> worker index in the arena.
        std::size_t thread_slot() const { return budget.inner == 1 ? 0 : housepi::parallel::current_thread_index(); }
        /// Return stream of chain `c`.
        /// * default -> chain 0, the one evaluation and simulation draw from.
        housepi::RndGen &prng(Chain c = 0) const { return prng_pool_[c].gen; }
        /// Restart the calling thread's household generator of chain `c`, return it.
        /// Its state is keyed_state() of (seed, `c`, `step`, `H`), whichever thread calls.
        housepi::RndGen &start_house_prng(Chain c, Step step, std::size_t H) const
        {
            housepi::RndGen &gen = house_prng(c);
            gen.restart(housepi::keyed_state(static_cast<std::uint64_t>(seed), c, step, H));
            return gen;
        }
        /// Return the calling thread's household generator of chain `c`.
        /// @note Draws continue the household of the last start_house_prng() on this thread.
        housepi::RndGen &house_prng(Chain c) const
        {
            return house_pool_[c % budget.outer + thread_slot() * budget.outer].gen;
        }

        /// Sample an exponential-distributed variable with rate 1, from the generator of chain `c`.
        double exponential(Chain c) const { return housepi::standard_exponential(prng(c)); }

        //////  Progress  //////

        /// Count one step.
        /// * stop flag raised -> throws parallel::Cancelled.
        void count_step()
        {
            if (stop_flag && stop_flag->load(std::memory_order_relaxed))
                throw parallel::Cancelled{};
            n_done.fetch_add(1, std::memory_order_relaxed);
        }
        /// Return fraction of the run's steps counted, from 0 to 1.
        double fraction_done() const;

    private:
        /// Generator on a cache line of its own, so threads drawing from neighbouring streams share no line.
        struct alignas(128) LineRndGen
        {
            housepi::RndGen gen;
        };
        /// Chain streams, seeded once: stream `c` starts `c * 2^128` draws past `seed`.
        mutable std::vector<LineRndGen> prng_pool_ = {};
        /// Household generators, one per thread, restarted at each household.
        /// * chain `c` at thread slot `s` -> generator `c % budget.outer + s * budget.outer`.
        mutable std::vector<LineRndGen> house_pool_ = {};

        /// Associate non-default user settings.
        void helper_non_default_build_(const housepi::str_table &input);
    };
}
