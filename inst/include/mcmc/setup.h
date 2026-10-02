/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/random.h"

//////////////////////////////  Setup  //////////////////////////////

namespace housepi
{
    class ProgressBar;

    /// Split workers between the two levels of parallel dispatch.
    struct Budget
    {
        std::size_t outer = 1; /// Workers spread across chains.
        std::size_t inner = 1; /// Workers spread within one chain.
    };

    /// Configure an MCMC run.
    ///
    /// - holds user, advanced and derived options;
    /// - owns the random number generators and the progress bar.
    struct Setup
    {
        //////  User config  //////

        std::int32_t seed;                    /// Seed of the pseudo-random number generators.
        unsigned n_threads;                   /// Number of threads for parallel computations.
        unsigned long n_steps;                /// Total number of MCMC iterations per chain.
        unsigned long n_burn;                 /// Number of initial burn-in iterations to discard.
        unsigned long n_thin;                 /// Thinning interval, storing every n_thin-th iteration.
        unsigned long n_adapt;                /// Number of adaptive MCMC steps, 0 for no adaptation.
        unsigned long n_chains;               /// Number of independent MCMC chains.
        unsigned long progress;               /// Progress bar update interval in seconds.
        std::vector<std::string> id_draws;    /// Individual IDs whose latent draws are traced.
        std::vector<std::string> house_draws; /// Household IDs whose infected members are traced.

        //////  Advanced config  //////

        unsigned long n_adapt_window; /// Effective window for the smoothed accept-rate.
        double time_step;             /// Real-time integration step for simulation.

        //////  Derived config  //////

        Budget budget;               /// Workers allotted to the chain and the per-chain dispatch level.
        unsigned long c_sample_size; /// Sample size for a single chain after thinning.
        double adapt_smoothing;      /// Smoothing factor of the accept-rate average, 1 / n_adapt_window.

        //////  Constructor  //////

        Setup();
        ~Setup();

        /// Configure progress bar for R sessions.
        /// @note Call once after construction when running from R.
        void init_R_context(unsigned barwidth, bool is_interactive,
                            std::function<void()> interrupt_func);

        //////  Setter  //////

        /// Set the number of threads, capped at those available.
        void set_n_threads();
        /// Reset every field to its default.
        void reset_defaults();
        /// Build setup from input table.
        void build(const housepi::str_table &input);
        /// Assert consistency of MCMC options and derive the sample size.
        void init_mcmc_configuration(bool is_adaptive);

        //////  Flags  //////

        /// Whether the run takes more than one step.
        bool is_mcmc_sample() const noexcept { return n_steps > 1; }

        //////  Random  //////

        /// Seed one generator per thread.
        ///
        /// - thread 0 takes `seed`;
        /// - each next thread jumps 2^128 steps ahead of the previous one.
        void set_prng_seed();
        /// Return generator of the calling thread.
        housepi::rnd_gen &prng() const;

        /// Sample a uniform-distributed variable in [0, 1).
        double uniform() const { return housepi::standard_uniform(prng()); }
        /// Sample an exponential-distributed variable with rate 1.
        double exponential() const { return housepi::standard_exponential(prng()); }

        //////  Progress bar  //////

        /// Update progress bar display text.
        void prog_set_text(std::string text = "");
        /// Increment progress bar by one iteration.
        void prog_update();
        /// Close the progress bar, returning the seconds it timed.
        double prog_close();
        /// Close progress bar on error `e` and raise it.
        [[noreturn]] void crash(const std::exception &e);

        /// Open progress bar for standard MCMC sample.
        void prog_open_for_mcmc_sample();
        /// Open progress bar for simulation.
        void prog_open_for_simulation();

    private:
        mutable std::vector<housepi::rnd_gen> prng_parallel_ = {}; /// One generator per thread.
        std::unique_ptr<housepi::ProgressBar> progressbar_;        /// Progress bar UI component.

        /// Associate non-default user setup.
        void helper_non_default_build_(const housepi::str_table &input);
    };
}
