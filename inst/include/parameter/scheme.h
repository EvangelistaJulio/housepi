/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <vector>

#include "tools/keys.h"
#include "parameter/definition.h"

//////////////////////////////  param::Scheme  //////////////////////////////

namespace housepi::param
{
    /// Sampling scheme of parameters: evaluation order, what a sweep samples and adapts.
    ///
    /// - `constant` + `sample`: every parameter;
    /// - `func`: subset of `constant`, deterministic function parameters;
    /// - `sample_standard` + `sample_cached`: `sample`.
    struct Scheme
    {
        std::vector<size_t> order = {};    /// Definition rows, producer before consumer.
        std::vector<size_t> constant = {}; /// Constant parameters.
        std::vector<size_t> func = {};     /// Function parameters.
        std::vector<size_t> sample = {};   /// Parameters to sample.
        std::vector<size_t> adaptive = {}; /// Parameters with an adaptive proposal, a subset of `sample`.

        std::vector<size_t> sample_standard = {}; /// Sampled parameters feeding no transmission cache.
        std::vector<size_t> sample_cached = {};   /// Sampled parameters feeding the transmission cache.

        std::vector<std::vector<double>> simulation_values = {}; /// Pinned values per repetition, per parameter.

        /// Classify definitions `def`, visited in `rows` order, over `keys`.
        void build(const DefTable &def,
                   const KeyHolder &keys,
                   std::vector<size_t> rows);
        /// Pin parameter `p` to `values`, one per simulation repetition, and drop it from sampling.
        void pin(size_t p, std::vector<double> values);
        /// Mark parameter `p` as a producer of the transmission cache.
        void mark_produces_transm_cache(size_t p);
    };
}
