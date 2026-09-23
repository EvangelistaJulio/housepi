/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "progressbar.h"
#include "setup.h"

namespace housepi
{
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** Setup : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Setup::Setup()
        : progressbar_(std::make_unique<ProgressBar>()) { reset_defaults(); }
    Setup::~Setup() = default;

    void Setup::init_R_context(
        unsigned barwidth,
        bool is_interactive,
        std::function<void()> interrupt_func)
    {
        progressbar_ = std::make_unique<ProgressBar>(
            barwidth,
            is_interactive,
            std::move(interrupt_func));
    }

    void Setup::reset_defaults()
    {
        // user configuration
        seed = 1;
        n_threads = 1;
        n_steps = 1;
        n_burn = 0;
        n_thin = 1;
        n_adapt = 0;
        n_chains = 1;
        progress = 1;
        id_draws.clear();
        house_draws.clear();
        // advanced configuration
        n_adapt_window = 50;
        time_step = 0.05;
        // derived configuration
        budget = {};
        c_sample_size = 0;
    }

    void Setup::build(const str_table &input)
    {
        reset_defaults();                 // single source of truth for defaults
        helper_non_default_build_(input); // read user input
        set_n_threads();                  // define number of threads
        set_prng_seed();                  // pseudo-random number generator seeding
    }

    void Setup::helper_non_default_build_(const housepi::str_table &input)
    {
        // check entries
        const std::vector<std::string> &headers = read::extract_keys(input);

        // Random Initialization of Parameters
        if (input.count("progress") && !na::is(input.at("progress")[0]))
            progress = str::toul(input.at("progress")[0]);

        // N. of Threads
        if (input.count("threads") && !na::is(input.at("threads")[0]))
            n_threads = str::toi(input.at("threads")[0], "threads");

        // PRNG Seed (enforce 32bit integer)
        if (input.count("seed") && !na::is(input.at("seed")[0]))
            seed = static_cast<std::int32_t>(str::toi(input.at("seed")[0], "seed"));

        // N. of Chains
        if (input.count("chains") && !na::is(input.at("chains")[0]))
            n_chains = str::toul(input.at("chains")[0], "chains");

        // N. of Steps
        if (input.count("steps") && !na::is(input.at("steps")[0]))
            n_steps = str::toul(input.at("steps")[0], "steps");

        // N. of Adaptive Steps
        if (input.count("adapt") && !na::is(input.at("adapt")[0]))
            n_adapt = str::toul(input.at("adapt")[0], "adpat");

        // EMA window for the smoothed acceptance-rate estimate.
        if (input.count("adapt_window") && !na::is(input.at("adapt_window")[0]))
            n_adapt_window = str::toul(input.at("adapt_window")[0], "adapt_window");

        // N. of Thinning Steps
        if (input.count("thin") && !na::is(input.at("thin")[0]))
            n_thin = str::toul(input.at("thin")[0], "thin");

        // N. of Steps to "Burn-In"
        if (input.count("burn") && !na::is(input.at("burn")[0]))
            n_burn = str::toul(input.at("burn")[0], "burn");

        // Individual ID for Latent Trace
        if (input.count("id_draws") && !input.at("id_draws")[0].empty())
            id_draws = str::unique(read::parse_list(input.at("id_draws")[0]));

        // Household ID for Latent Trace
        if (input.count("house_draws") && !input.at("house_draws")[0].empty())
            house_draws = str::unique(read::parse_list(input.at("house_draws")[0]));

        // Real-time integration step for simulation
        if (input.count("time_step") && !na::is(input.at("time_step")[0]))
        {
            time_step = str::tod(input.at("time_step")[0], "time_step");
            // UPDATE: by symmetry it should be tested where appropriate
            if (time_step <= 0.)
                msg::stop("Setup 'time_step' must be strictly positive.");
        }
    }

    void Setup::set_n_threads()
    {
        n_threads = std::max(1u, n_threads);

#ifdef _STANDALONE
        if (n_threads != 1)
            msg::warn(
                "Multi-threading is disabled.\nRcppParallel is not available in standalone mode.");
        n_threads = 1;

#else
        size_t max_threads = housepi::parallel::max_num_threads();
        if (n_threads > max_threads)
        {
            msg::warn(
                "Reducing " + msg::arg("threads") + " to " +
                std::to_string(max_threads) + " available threads.");
            n_threads = max_threads;
        }
        housepi::parallel::set_num_threads(n_threads);
#endif
        // one dispatch level: every thread goes to the households of a single chain
        budget = {1, n_threads};
    }

    /// Assign seed for random generators.
    /// Concomitantely create different generators for different parallel threads.
    void Setup::set_prng_seed()
    {
        seed = std::max(1, seed); // NOTE: default seed for standalone
        prng_parallel_.resize(n_threads);
        prng_parallel_[0].seed(seed);
        for (size_t t = 1; t < n_threads; ++t)
        {
            prng_parallel_[t] = prng_parallel_[t - 1];
            prng_parallel_[t].jump();
        }
    }

    /// Return pseudo random number generator.
    rnd_gen &Setup::prng() const
    {
        return prng_parallel_[housepi::parallel::current_thread_index()];
    }

    void Setup::init_mcmc_configuration(bool is_adaptive)
    {
        // clip options mistakenly set to 0 has 1
        n_steps = std::max(1ul, n_steps);
        n_chains = std::max(1ul, n_chains);
        n_thin = std::max(1ul, n_thin);
        n_adapt_window = std::max(1ul, n_adapt_window);

        if (n_burn >= n_steps)
            msg::stop(
                "Wrong setup: " + msg::arg("burn") + " must be smaller than " + msg::arg("steps") + ".");
        if (n_thin > n_steps)
            msg::stop(
                "Wrong setup: " + msg::arg("thin") + "cannot be greater than " + msg::arg("steps") + ".");

        c_sample_size = (n_steps - n_burn) / n_thin;
        // adjust steps and sample size
        c_sample_size += (n_steps - n_burn) % n_thin > 0;
        n_steps = n_burn + c_sample_size * n_thin;

        // drop adaptation steps if nothing to adapt
        // if burn = 0, repurpose adaptation budget as burn
        if (!is_adaptive)
        {
            if (n_adapt > 0 && n_burn == 0)
            {
                n_burn = n_adapt;
                n_steps += n_adapt;
            }
            n_adapt = 0;
        }
    }

    void Setup::prog_set_text(std::string text)
    {
        if (na::is(text))
            text = "sample";
        progressbar_->set_text(text);
    }
    void Setup::prog_open_for_mcmc_sample()
    {
        prog_set_text();
        progressbar_->open((n_adapt + n_steps) * n_chains, progress);
    }
    void Setup::prog_open_for_simulation()
    {
        prog_set_text("simulate");
        progressbar_->open(n_steps, progress);
    }

    void Setup::prog_update() { progressbar_->update(); }
    double Setup::prog_close() { return progressbar_->close(); }
    [[noreturn]] void Setup::crash(const std::exception &e)
    {
        progressbar_->crash();
        msg::stop(e.what());
    }
}
