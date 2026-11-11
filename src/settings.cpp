/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "tools/string.h"
#include "mcmc/settings.h"

namespace housepi
{
    //////////////////////////////  Settings  //////////////////////////////

    void Settings::reset_defaults()
    {
        // user configuration
        seed = 1;
        n_threads = 1;
        n_steps = 1;
        n_thin = 1;
        n_adapt = 0;
        n_chains = 1;
        id_draws.clear();
        house_draws.clear();
        // advanced configuration
        time_step = 0.05;
        // derived configuration
        budget = {};
        c_sample_size = 0;
        n_burn = 0;
    }

    void Settings::build(const str_table &input)
    {
        reset_defaults();
        helper_non_default_build_(input);
        set_n_threads();
        set_prng_seed();
    }

    void Settings::helper_non_default_build_(const housepi::str_table &input)
    {
        // number of threads
        if (input.count("threads") && !na::is(input.at("threads")[0]))
            n_threads = str::toi(input.at("threads")[0], "threads");

        // seed, cast to a 32-bit integer
        if (input.count("seed") && !na::is(input.at("seed")[0]))
            seed = static_cast<std::int32_t>(str::toi(input.at("seed")[0], "seed"));

        // number of chains
        if (input.count("chains") && !na::is(input.at("chains")[0]))
            n_chains = str::toull(input.at("chains")[0], "chains");

        // number of steps
        if (input.count("steps") && !na::is(input.at("steps")[0]))
            n_steps = str::toull(input.at("steps")[0], "steps");

        // number of adaptive steps
        if (input.count("adapt") && !na::is(input.at("adapt")[0]))
            n_adapt = str::toull(input.at("adapt")[0], "adapt");

        // thinning interval
        if (input.count("thin") && !na::is(input.at("thin")[0]))
            n_thin = str::toull(input.at("thin")[0], "thin");

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
                msg::stop("Settings 'time_step' must be strictly positive.");
        }
    }

    void Settings::set_n_threads()
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
        // every thread to households, split across chains by init_mcmc_configuration()
        budget = {1, n_threads};
    }

    void Settings::set_prng_seed()
    {
        seed = std::max(1, seed);
        prng_pool_.clear();
        RndGen next(static_cast<RndGen::result_type>(seed));
        for (size_t k = 0; k < n_chains * n_threads; ++k)
        {
            prng_pool_.push_back({next});
            next.jump();
        }
    }

    void Settings::init_mcmc_configuration(bool is_adaptive)
    {
        // every state runs as many chains on as many threads, none idle
        // @note threads not a multiple of chains go to households
        // where a model may scale worse than over chains
        budget.outer = std::gcd<size_t>(n_chains, n_threads);
        budget.inner = n_threads / budget.outer;

        c_sample_size = n_steps / n_thin;

        // nothing to adapt: adaptation steps run as burn-in
        if (!is_adaptive)
        {
            n_burn = n_adapt;
            n_adapt = 0;
        }
    }

    double Settings::fraction_done() const
    {
        return static_cast<double>(n_done.load(std::memory_order_relaxed)) / ((n_adapt + n_burn + n_steps) * n_chains);
    }
}
