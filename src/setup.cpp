/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "tools/string.h"
#include "mcmc/setup.h"

namespace housepi
{
    //////////////////////////////  Setup  //////////////////////////////

    Setup::Setup() { reset_defaults(); }

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
        id_draws.clear();
        house_draws.clear();
        // advanced configuration
        n_adapt_window = 50;
        time_step = 0.05;
        // derived configuration
        budget = {};
        c_sample_size = 0;
        adapt_smoothing = 1.;
    }

    void Setup::build(const str_table &input)
    {
        reset_defaults();                 // defaults
        helper_non_default_build_(input); // user input
        set_n_threads();                  // number of threads
        set_prng_seed();                  // generator seeds
    }

    void Setup::helper_non_default_build_(const housepi::str_table &input)
    {
        // number of threads
        if (input.count("threads") && !na::is(input.at("threads")[0]))
            n_threads = str::toi(input.at("threads")[0], "threads");

        // seed, cast to a 32-bit integer
        if (input.count("seed") && !na::is(input.at("seed")[0]))
            seed = static_cast<std::int32_t>(str::toi(input.at("seed")[0], "seed"));

        // number of chains
        if (input.count("chains") && !na::is(input.at("chains")[0]))
            n_chains = str::toul(input.at("chains")[0], "chains");

        // number of steps
        if (input.count("steps") && !na::is(input.at("steps")[0]))
            n_steps = str::toul(input.at("steps")[0], "steps");

        // number of adaptive steps
        if (input.count("adapt") && !na::is(input.at("adapt")[0]))
            n_adapt = str::toul(input.at("adapt")[0], "adapt");

        // window of the smoothed acceptance rate
        if (input.count("adapt_window") && !na::is(input.at("adapt_window")[0]))
            n_adapt_window = str::toul(input.at("adapt_window")[0], "adapt_window");

        // thinning interval
        if (input.count("thin") && !na::is(input.at("thin")[0]))
            n_thin = str::toul(input.at("thin")[0], "thin");

        // burn-in steps
        if (input.count("burn") && !na::is(input.at("burn")[0]))
            n_burn = str::toul(input.at("burn")[0], "burn");

        // individual IDs whose latent draws are traced
        if (input.count("id_draws") && !input.at("id_draws")[0].empty())
            id_draws = str::unique(str::split_list(input.at("id_draws")[0]));

        // household IDs whose infected members are traced
        if (input.count("house_draws") && !input.at("house_draws")[0].empty())
            house_draws = str::unique(str::split_list(input.at("house_draws")[0]));

        // integration step of simulation
        if (input.count("time_step") && !na::is(input.at("time_step")[0]))
        {
            time_step = str::tod(input.at("time_step")[0], "time_step");
            if (time_step <= 0.)
                msg::stop("Setup 'time_step' must be strictly positive.");
        }
    }

    void Setup::set_n_threads()
    {
        n_threads = std::max(1u, n_threads);

        size_t max_threads = housepi::parallel::max_num_threads();
        if (n_threads > max_threads)
        {
            msg::warn(
                "Reducing " + msg::arg("threads") + " to " +
                std::to_string(max_threads) + " available threads.");
            n_threads = max_threads;
        }
        // every thread goes to the households of a single chain
        budget = {1, n_threads};
    }

    void Setup::set_prng_seed()
    {
        seed = std::max(1, seed); // seeds below 1 fall back to 1
        prng_parallel_.resize(n_threads);
        prng_parallel_[0].seed(seed);
        for (size_t t = 1; t < n_threads; ++t)
        {
            prng_parallel_[t] = prng_parallel_[t - 1];
            prng_parallel_[t].jump();
        }
    }

    rnd_gen &Setup::prng() const
    {
        return prng_parallel_[budget.inner == 1 ? 0 : housepi::parallel::current_thread_index()];
    }

    void Setup::init_mcmc_configuration(bool is_adaptive)
    {
        adapt_smoothing = 1. / n_adapt_window;
        c_sample_size = n_steps / n_thin;

        // if nothing to adapt, adaptation steps are:
        // - converted in burn-in steps when burn is 0;
        // - dropped otherwise.
        if (!is_adaptive)
        {
            if (n_burn == 0)
                n_burn = n_adapt;
            n_adapt = 0;
        }
    }

    void Setup::count_iter()
    {
        if (stop_flag && stop_flag->load(std::memory_order_relaxed))
            throw parallel::Cancelled{};
        n_iter.store(n_iter.load(std::memory_order_relaxed) + 1, std::memory_order_relaxed);
    }
    double Setup::fraction_iter() const
    {
        return static_cast<double>(n_iter.load(std::memory_order_relaxed)) / ((n_adapt + n_burn + n_steps) * n_chains);
    }
}
