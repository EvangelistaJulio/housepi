/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>
#include <vector>

#include "tools/dictionary.h"
#include "tools/random.h"
#include "data/individual.h"
#include "data/status.h"

namespace housepi::data
{
    struct Wide;
}

//////////////////////////////  StatusProcess  //////////////////////////////

namespace housepi
{
    /// Infected status drawn per member, weighted by a rule table over covariates and status.
    ///
    /// - weight of infected level `k`: product over rule groups of the most specific matching rule;
    /// - members of one covariate class share one block.
    struct StatusProcess
    {
        /// Resolve the weights of every member from rule table `input`, over `n_levels` infected levels.
        /// - empty `input`: no status process.
        /// @note Coverage, ties and positive totals are asserted R-side.
        void build(const data::Wide &data_wide,
                   const str_table &input,
                   size_t n_levels);

        /// Whether a status process is defined.
        bool is_defined() const noexcept { return !blocks_.empty(); }

        /// Draw the infected status of member `indiv`.
        /// - no status process: the lowest infected status.
        StatusLvl draw(const Individual &indiv, rnd_gen &prng) const;

        /// Infected statuses member `indiv` can draw, those of positive weight.
        /// - no status process: the lowest infected status.
        std::vector<StatusLvl> drawable(const Individual &indiv) const;

    private:
        std::vector<std::vector<double>> blocks_ = {}; /// Weights per covariate class, one per infected level.
        std::vector<double> totals_ = {};              /// Sum of the weights of every block.
        std::vector<size_t> block_of_ = {};            /// Block of every member, its covariate class.
    };
}
