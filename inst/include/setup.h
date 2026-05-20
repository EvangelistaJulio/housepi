/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <functional>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "tools/random.h"
#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Setup : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    class ProgressBar;

    /// Configuration and control parameters for MCMC sampling.
    /// Manages random number generation and progress tracking.
    struct Setup
    {
        ///////////////////////////
        //////  USER CONFIG  //////
        ///////////////////////////

        /// Expected header keys for configuration input file.
        const std::vector<std::string> file_headers = {"steps", "burn", "thin",
                                                       "adapt", "chains", "threads",
                                                       "id_draws", "seed", "progress"};

        std::int32_t seed = 1;      /// Seed for pseudo-random number generator initialization.
        unsigned n_threads = 1;     /// Number of OpenMP threads for parallelization.
        unsigned long n_steps = 1;  /// Total number of MCMC iterations per chain.
        unsigned long n_burn = 0;   /// Number of initial burn-in iterations to discard.
        unsigned long n_thin = 1;   /// Thinning interval (store every n_thin-th iteration).
        unsigned long n_adapt = 0;  /// Number of adaptive MCMC steps (0 = no adaptation).
        unsigned long n_chains = 1; /// Number of independent MCMC chains to run in parallel.
        unsigned long progress = 1; /// Progress bar update interval in seconds.
        unsigned long n_rungs = 0;  /// Number of rungs for power posterior analysis (0 = standard MCMC).

        //////////////////////////////
        //////  DERIVED CONFIG  //////
        //////////////////////////////

        unsigned long chain_id = 1;             /// Current active MCMC chain (1-indexed).
        unsigned long rung_id = 0;              /// Current active rung in power posterior analysis (0-indexed).
        unsigned long c_sample_size = 0;        /// Sample size for a single chain after thinning.
        unsigned long c_n_sweeps = 0;           /// Total number of sweeps per chain.
        double inv_T = 1.;                      /// Inverse temperature for power posterior (1.0 = posterior).
        std::vector<std::string> id_draws = {}; /// IDs of augmented variables to trace and save.

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        Setup();
        ~Setup();

        /// Configure progress bar for R sessions.
        /// Must be called once after construction when running from R.
        void init_R_context(unsigned barwidth, bool is_interactive,
                            std::function<void()> interrupt_func);

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Set the number of threads for parrell computations (default is 1).
        void set_n_threads();

        /// Build Setup configuration from input table.
        void build(const housepi::str_table &input);

        /// Validate MCMC setup parameters for consistency and feasibility.
        void assert_mcmc_setup();
        /// Remove adaptive phase parameters after burn-in completion.
        void drop_mcmc_adapt();

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Check if running standard sample.
        bool is_mcmc_sample() const noexcept { return n_steps > 1; }
        /// Check if running power posterior analysis.
        bool is_power_post() const noexcept { return n_rungs > 0; }

        //////////////////////
        //////  RANDOM  //////
        //////////////////////

        /// Initialize random number generator with configured seed.
        void set_prng_seed();
        /// Get reference to the active pseudo-random number generator.
        housepi::rnd_gen &prng() const;

        /// Sample a uniform-distributed variable in [0, 1).
        double uniform() const { return uniform_(prng()); }
        /// Sample an exponential-distributed variable with rate 1.
        double exponential() const { return exponential_(prng()); }

        ////////////////////////////
        //////  PROGRESS-BAR  //////
        ////////////////////////////

        /// Update progress bar display text.
        void prog_set_text(std::string text = "");

        /// Open progress bar for standard MCMC sample.
        void prog_open_for_mcmc_sample();
        /* --FUTURE-LOG-EVIDENCE--
        /// Open progress bar for power posterior analysis.
        void prog_open_for_power_post();
        */
        /// Open progress bar for simulation.
        void prog_open_for_simulation();

        /// Increment progress bar by one iteration.
        void prog_update();
        /// Close and finalize progress bar display.
        void prog_close();
        /// Mark progress bar as crashed/interrupted.
        [[noreturn]] void crash(const std::exception &e);

    private:
        mutable std::vector<housepi::rnd_gen> prng_parallel_ = {};  /// Thread-local PRNGs for parallel chains.
        mutable std::uniform_real_distribution<double> uniform_;    /// Uniform [0, 1) distribution sampler.
        mutable std::exponential_distribution<double> exponential_; /// Exponential distribution sampler.
        std::unique_ptr<housepi::ProgressBar> progressbar_;         /// Progress bar UI component.

        /// Associate non-default user setup.
        void helper_non_default_build_(const housepi::str_table &input);
    };
}
