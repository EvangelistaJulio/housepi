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
    /// * `constant` + `sample`: every parameter;
    /// * `func`: subset of `constant`, deterministic function parameters;
    /// * `sample_standard` + `sample_cached`: `sample`.
    struct Scheme
    {
        /// Definition rows, producer before consumer.
        std::vector<size_t> order = {};
        std::vector<size_t> constant = {};
        std::vector<size_t> func = {};
        std::vector<size_t> sample = {};
        /// Parameters with an adaptive proposal, a subset of `sample`.
        std::vector<size_t> adaptive = {};

        /// Sampled parameters feeding no transmission cache.
        std::vector<size_t> sample_standard = {};
        /// Sampled parameters feeding the transmission cache.
        std::vector<size_t> sample_cached = {};

        /// Pinned values per repetition, per parameter.
        /// * not pinned -> null.
        std::vector<const std::vector<double> *> simulation_values = {};

        /// Classify definitions `def`, visited in `rows` order, over `keys`.
        void build(const DefTable &def,
                   const KeyHolder &keys,
                   std::vector<size_t> rows);
        /// Pin parameter `p` to `values`, one per simulation repetition, and drop it from sampling.
        /// @note The pin points into `values`, which outlives it.
        void pin(size_t p, const std::vector<double> &values);
        /// Mark parameter `p` as a producer of the transmission cache.
        void mark_produces_transm_cache(size_t p);
    };
}
