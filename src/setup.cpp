/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <chrono>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/progress.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "setup.h"

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** Setup : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

Setup::Setup()
    : uniform_(0., 1.),
      progressbar_(std::make_unique<hsp::progress::Bar>()) {}
Setup::~Setup() = default;

void Setup::build(const hsp::str_table &input)
{
    const auto time_raw = std::chrono::steady_clock::now().time_since_epoch().count();
    // ensure time_now is 64bit
    const auto time_now = static_cast<std::uint64_t>(time_raw);
    // xor between upper and lower part of time stamp
    const std::uint64_t mixed = time_now ^ (time_now >> 32);
    // default seed for pseudo random number generator
    // note: constrain to 32bit for compatibility with R (mask sign bit)
    seed = static_cast<std::int32_t>(mixed & 0x7FFFFFFF);

    chain = 1;        // current chain
    rung = 0;         // current thermodynamic rung
    id_draws.clear(); // drop previous id

    if (!input.empty())
    {
        // check no multiple entries if internal error
        if (hsp::read::csv_n_rows(input) != 1)
            throw std::runtime_error(
                "Internal error: multiple entries found in MCMC setup definition.");

        // check entries
        const std::vector<std::string> &headers = hsp::read::extract_keys(input);
        hsp::str::assert_no_unused_entries(headers, file_headers);

        // Random Initialization of Parameters
        if (input.count("progress") && !hsp::na::is(input.at("progress")[0]))
            progress = hsp::str::toul(input.at("progress")[0]);

        // N. of Threads
        if (input.count("threads") && !hsp::na::is(input.at("threads")[0]))
            n_threads = hsp::str::toi(input.at("threads")[0], "threads");

        // PRNG Seed
        if (input.count("seed") && !hsp::na::is(input.at("seed")[0]))
            seed = hsp::str::toi32(input.at("seed")[0], "seed");

        // N. of Chains
        if (input.count("chains") && !hsp::na::is(input.at("chains")[0]))
            n_chains = hsp::str::toul(input.at("chains")[0], "chains");

        // N. of Steps
        if (input.count("steps") && !hsp::na::is(input.at("steps")[0]))
            n_steps = hsp::str::toul(input.at("steps")[0], "steps");

        // N. of Adaptive Steps
        if (input.count("adapt") && !hsp::na::is(input.at("adapt")[0]))
            n_adapt = hsp::str::toul(input.at("adapt")[0], "adpat");

        // N. of Thinning Steps
        if (input.count("thin") && !hsp::na::is(input.at("thin")[0]))
            n_thin = hsp::str::toul(input.at("thin")[0], "thin");

        // N. of Steps to "Burn-In"
        if (input.count("burnin") && !hsp::na::is(input.at("burnin")[0]))
            n_burnin = hsp::str::toul(input.at("burnin")[0], "burnin");

        // Individual ID for Augmented Trace
        if (input.count("id_draws") && !hsp::na::is(input.at("id_draws")[0]))
            id_draws = hsp::read::parse_list(input.at("id_draws")[0]);
    }

    set_n_threads(); // define number of threads
    set_prng_seed(); // pseudo-random number generator seeding
}

/// Set the number of threads for parrell computations (default is 1).
void Setup::set_n_threads()
{
    if (n_threads == 0)
        throw std::invalid_argument(
            "N. of threads was not greater than 0.");

#ifdef _OPENMP
    unsigned int max_threads = omp_get_num_procs();

    if (n_threads > max_threads)
    {
        hsp::msg::warning(
            "N. of threads reduced to " + std::to_string(max_threads) + " available threads.");
        n_threads = max_threads;
    }

    // declare the global number of omp threads
    omp_set_num_threads(n_threads);

#else
    if (n_threads != 1)
        hsp::msg::warning(
            "Multi-threading was disabled. OpenMP is not supported on current machine.");

    n_threads = 1;
#endif
}

/// Assign seed for random generators.
/// Concomitantely create different generators for different parallel threads
void Setup::set_prng_seed()
{
    hsp::rnd_gen aux_prng;
    aux_prng.seed(seed);
    std::uniform_int_distribution<unsigned long> seed_generator(0, 0XFFFFFFFF);

    prng_parallel_.resize(n_threads);
    for (unsigned n = 0; n < n_threads; ++n)
        prng_parallel_[n].seed(seed_generator(aux_prng));
}

/// Return pseudo random number generator.
hsp::rnd_gen &Setup::prng()
{
#ifdef _OPENMP
    return prng_parallel_[omp_get_thread_num()];
#else
    return prng_parallel_[0];
#endif
}

void Setup::assert_mcmc_setup()
{
    if (inv_T < 0. || inv_T > 1.)
        throw std::range_error(
            "Inverse temperature " + std::to_string(inv_T) + " was out of range [0, 1].");

    if (n_chains < 1)
        throw std::out_of_range(
            "Number of chains was not greater than 0.");

    // Interpret thinning equal to 0 has been interpreted as 1.
    if (n_thin == 0)
        n_thin = 1;

    if (n_steps < 1)
        throw std::out_of_range(
            "N. of steps was not greater than 0.");
    if (n_burnin >= n_steps)
        throw std::out_of_range(
            "Burn-in must was not smaller than n. of steps.");
    if (n_thin > n_steps)
        throw std::out_of_range(
            "Thinning was not be smaller than (or equal to) n. of steps. ");

    c_sample_size = (n_steps - n_burnin) / n_thin;
    // adjust steps and sample size
    c_sample_size += (n_steps - n_burnin) % n_thin > 0;
    n_steps = n_burnin + c_sample_size * n_thin;

    if (c_sample_size == 0)
        throw std::out_of_range(
            "Choice of steps, burn-in and thin generated an empty sample.");

    // assign number of sweeps
    // NOTE: because of rounding I can't put all in the same line
    c_n_sweeps = n_adapt / n_thin;
    c_n_sweeps += n_burnin / n_thin;
    c_n_sweeps += c_sample_size;
}

// Drop adaptation steps.
void Setup::drop_mcmc_adapt()
{
    if (n_adapt > 0)
        hsp::msg::warning(
            "Dropped adaptation steps: no parametric MCMC proposals were found.");
    n_adapt = 0;

    // re-assign number of sweeps
    c_n_sweeps = n_burnin / n_thin;
    c_n_sweeps += c_sample_size;
}

/// The text to print in the Markov chian progress bar.
std::string Setup::fmt_progtxt()
{
    if (is_power_post())
        return "power post " + std::to_string(rung) + "/" + std::to_string(n_rungs);
    if (n_chains > 1)
        return "sample " + std::to_string(chain) + "/" + std::to_string(n_chains);
    return "sample";
}

void Setup::prog_declare(const std::string &text)
{
    prog_set_text(text);
    progressbar_->set_thintime(progress);
}
void Setup::prog_set_text(const std::string &text) { progressbar_->set_text(hsp::na::is(text) ? fmt_progtxt() : text); }
void Setup::prog_open_for_inference() { progressbar_->open(c_n_sweeps * n_chains); }
void Setup::prog_open_for_power_post() { progressbar_->open(c_n_sweeps * n_chains * (n_rungs + 1)); }
void Setup::prog_update() { progressbar_->update(); }
void Setup::prog_close() { progressbar_->close(); }
void Setup::prog_crash() { progressbar_->crash(); }
