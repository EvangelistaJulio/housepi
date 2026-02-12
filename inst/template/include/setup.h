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
}

struct Setup
{
    Setup();
    ~Setup();

    const std::vector<std::string> file_headers = {"steps", "burnin", "thin",
                                                   "adapt", "chains", "threads",
                                                   "id_draws", "seed", "progress"};

    unsigned long n_steps = 1;  /// n. iterations
    unsigned long n_burnin = 0; /// n. burnin iterations
    unsigned long n_thin = 1;   /// n. thinning iterations
    unsigned long n_adapt = 0;  /// n. (finite) adaptive steps
    unsigned long n_chains = 1; /// n. chains
    unsigned n_threads = 1;     /// n. of OMP threads
    std::int32_t seed;          /// random generator seed
    unsigned long progress = 1; /// thinning seconds for progress bar
    unsigned long n_rungs = 0;  /// n. of rungs

    unsigned long c_sample_size = 0;        /// single chain sample size
    unsigned long c_n_sweeps = 0;           /// single chain n. of sweeps
    double inv_T = 1.;                      /// inverse temperature
    std::vector<std::string> id_draws = {}; /// augmented trace to save

    unsigned long chain = 1; /// current MCMC chain
    unsigned long rung = 0;  /// current MCMC thermodynamic rung

    bool is_power_post() const { return n_rungs > 0; }
    std::string fmt_progtxt();

    void set_n_threads();
    void build(const hsp::str_table &input);
    void assert_mcmc_setup();
    void drop_mcmc_adapt();

    // progress bar
    void prog_declare(const std::string &text = "");
    void prog_set_text(const std::string &text = "");
    void prog_open_for_inference();
    void prog_open_for_power_post();
    void prog_update();
    void prog_close();
    void prog_crash();

    /// Set seed for pseudo random number generator.
    void set_prng_seed();
    /// Return pseudo random number generator.
    hsp::rnd_gen &prng();
    /// Sample a uniform-distributed variable between 0 and 1.
    double uniform() { return uniform_(prng()); }
    /// Sample a exponential-distributed variable between 0 and Inf.
    double exponential() { return exponential_(prng()); }

private:
    std::vector<hsp::rnd_gen> prng_parallel_ = {};      // pseudo-random number generators
    std::uniform_real_distribution<double> uniform_;    // uniform sampler
    std::exponential_distribution<double> exponential_; // exponential sampler
    std::unique_ptr<hsp::ProgressBar> progressbar_;     // progress bar
};
