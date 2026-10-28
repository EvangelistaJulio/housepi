/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>
#include <vector>

#include "tools/nan.h"

//////////////////////////////  rates::Anchor  //////////////////////////////

namespace housepi
{
    struct Parameter;
}

namespace housepi::rates
{
    /// Share of the units of its reference from which a rate moves with its baseline.
    constexpr double BLOCK_MIN_SHARE = 0.1;

    /// Baseline parameter and the rate ratios scaling it, as declared by the model.
    struct Anchor
    {
        size_t baseline;           /// Parameter scaled by the rates.
        std::vector<size_t> rates; /// Group-rate keys scaling the baseline.
        bool contacts;             /// Whether contact rates scale the same term.
    };

    /// Scale anchoring of one rule group, read off its rules (see Rates::assert_identifiability).
    ///
    /// The rates are kept, not counted: whether they move depends on the pinning of the run.
    /// @note Every member or pair reading the rate matches a rule of the group, asserted R-side.
    struct ScaleAnchor
    {
        bool anchored = false;                   /// Whether a rule of the group names the reference.
        std::vector<size_t> scaling_params = {}; /// Distinct rates the group applies.
        double n_reference = 0.;                 /// Units the reference rules resolve, NA if any is unknown.
        std::vector<double> n_params = {};       /// Units the rules of each of `scaling_params` resolve, NA if any is unknown.

        /// List the rates holding at least a minimum share of the units of the reference.
        /// * no reference, or unknown units -> skipped.
        std::vector<size_t> list_heavy() const
        {
            std::vector<size_t> out = {};
            if (!anchored || na::is(n_reference) || !(n_reference > 0.))
                return out;
            for (size_t i = 0; i < scaling_params.size(); ++i)
                if (!na::is(n_params[i]) && n_params[i] >= BLOCK_MIN_SHARE * n_reference)
                    out.push_back(scaling_params[i]);
            return out;
        }

        /// Whether the overall scale of the rule group is free.
        /// * a rule names the reference -> not free;
        /// * some rate applied is pinned -> not free, it fixes the level of the whole group.
        bool is_free(const Parameter &param) const;
    };

    /// Count the rule groups in `anchors` whose overall scale is free.
    size_t count_scale_free(const std::vector<ScaleAnchor> &anchors,
                            const Parameter &param);
}
