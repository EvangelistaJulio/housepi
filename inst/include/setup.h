/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

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

namespace hsp
{
    class ProgressBar;

    /// Configuration and control parameters for MCMC sampling.
    /// Manages random number generation and progress tracking.
    struct Setup
    {
        /// Expected header keys for configuration input file.
        const std::vector<std::string> file_headers = {"steps", "burnin", "thin",
                                                       "adapt", "chains", "threads",
                                                       "id_draws", "seed", "progress"};

        unsigned long n_steps = 1;  /// Total number of MCMC iterations per chain.
        unsigned long n_burnin = 0; /// Number of initial burn-in iterations to discard.
        unsigned long n_thin = 1;   /// Thinning interval (store every n_thin-th iteration).
        unsigned long n_adapt = 0;  /// Number of adaptive MCMC steps (0 = no adaptation).
        unsigned long n_chains = 1; /// Number of independent MCMC chains to run in parallel.
        unsigned n_threads = 1;     /// Number of OpenMP threads for parallelization.
        std::int32_t seed;          /// Seed for pseudo-random number generator initialization.
        unsigned long progress = 1; /// Progress bar update interval in seconds.
        unsigned long n_rungs = 0;  /// Number of rungs for power posterior analysis (0 = standard MCMC).

        unsigned long c_sample_size = 0;        /// Sample size for a single chain after thinning.
        unsigned long c_n_sweeps = 0;           /// Total number of sweeps per chain.
        double inv_T = 1.;                      /// Inverse temperature for power posterior (1.0 = posterior).
        std::vector<std::string> id_draws = {}; /// IDs of augmented variables to trace and save.

        unsigned long chain = 1; /// Current active MCMC chain (1-indexed).
        unsigned long rung = 0;  /// Current active rung in power posterior analysis (0-indexed).

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        Setup();
        ~Setup();

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Set the number of threads for parrell computations (default is 1).
        void set_n_threads();

        /// Build Setup configuration from input table.
        /// @param input Key-value pairs matching `file_headers`.
        void build(const hsp::str_table &input);

        /// Validate MCMC setup parameters for consistency and feasibility.
        /// @throws std::out_of_range if parameters are invalid.
        void assert_mcmc_setup();

        /// Remove adaptive phase parameters after burn-in completion.
        void drop_mcmc_adapt();

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Check if running power posterior analysis.
        /// @return true if n_rungs > 0 (multiple temperature rungs configured).
        bool is_power_post() const { return n_rungs > 0; }

        //////////////////////
        //////  RANDOM  //////
        //////////////////////

        /// Initialize random number generator with configured seed.
        void set_prng_seed();

        /// Get reference to the active pseudo-random number generator.
        /// @return Reference to PRNG for current thread.
        hsp::rnd_gen &prng() const;

        /// Sample a uniform-distributed variable in [0, 1).
        /// @return Uniform random value.
        double uniform() const { return uniform_(prng()); }

        /// Sample an exponential-distributed variable with rate 1.
        /// @return Exponential random value in [0, inf).
        double exponential() const { return exponential_(prng()); }

        ////////////////////////////
        //////  PROGRESS-BAR  //////
        ////////////////////////////

        /// Declare and initialize progress bar with optional label.
        /// @param text Optional description to display.
        void prog_declare(const std::string &text = "");

        /// Update progress bar display text.
        /// @param text New text to display.
        void prog_set_text(const std::string &text = "");

        /// Open progress bar for standard MCMC inference.
        void prog_open_for_inference();

        /// Open progress bar for power posterior analysis.
        void prog_open_for_power_post();

        /// Increment progress bar by one iteration.
        void prog_update();

        /// Close and finalize progress bar display.
        void prog_close();

        /// Mark progress bar as crashed/interrupted.
        void prog_crash();

        /// Formatted progress text for display.
        /// @return Formatted string with current MCMC/power posterior state.
        std::string fmt_progtxt();

    private:
        mutable std::vector<hsp::rnd_gen> prng_parallel_ = {};      /// Thread-local PRNGs for parallel chains.
        mutable std::uniform_real_distribution<double> uniform_;    /// Uniform [0, 1) distribution sampler.
        mutable std::exponential_distribution<double> exponential_; /// Exponential distribution sampler.
        std::unique_ptr<hsp::ProgressBar> progressbar_;             /// Progress bar UI component.

        /// Associate non-default user setup.
        void helper_non_default_build_(const hsp::str_table &input);
    };
}
