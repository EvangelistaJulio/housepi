/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <vector>

#include "tools/keys.h"
#include "tools/random.h"
#include "parameter/scheme.h"
#include "variable.h"

//////////////////////////////  param::State  //////////////////////////////

namespace housepi::param
{
    /// Mutable per-run storage of parameters, moved by inference and simulation.
    struct State
    {
        /// Parameters, by key.
        std::vector<Variable> vec = {};

        /// Allocate one Variable per key of `keys`, named after it.
        void allocate(const KeyHolder &keys);
    };

    //////  Sampling  //////

    /// Initialize parameters to sample, for inference.
    void restart(const Scheme &scheme,
                 State &state,
                 RndGen &prng);
    /// Draw parameters to sample from their prior, for simulation.
    /// Unlike restart(), it uses no proposal.
    inline void draw_sample(const Scheme &scheme,
                            State &state,
                            RndGen &prng)
    {
        for (size_t p : scheme.sample)
            state.vec[p].draw_value(prng);
    }
    /// Load pinned values of simulation repetition `r`.
    void load_sample_values(const Scheme &scheme,
                            State &state,
                            size_t r);
    /// Finalize every parameter (see Variable::finalize).
    void finalize(const Scheme &scheme,
                  State &state);
}
