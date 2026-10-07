/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
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
        std::vector<Variable> vec = {}; /// Parameters, by key.

        /// Allocate one Variable per key of `khold`, named after it.
        void allocate(const KeyHolder &khold);
    };

    //////  Sampling  //////

    /// Initialize parameters to sample, for inference.
    void restart(const Scheme &scheme,
                 State &state,
                 rnd_gen &prng);
    /// Draw parameters to sample from their prior, for simulation.
    /// Unlike restart(), it uses no proposal.
    inline void draw_sample(const Scheme &scheme,
                            State &state,
                            rnd_gen &prng)
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
