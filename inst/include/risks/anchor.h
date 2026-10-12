/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>
#include <vector>

//////////////////////////////  risks::Anchor  //////////////////////////////

namespace housepi
{
    struct Parameter;
}

namespace housepi::risks
{
    /// Baseline parameter and the relative risks scaling it, as declared by the model.
    struct Anchor
    {
        size_t baseline;           /// Parameter scaled by the risks.
        std::vector<size_t> risks; /// Group-risk keys scaling the baseline.
        bool contacts;             /// Whether contact risks scale the same term.
    };

    /// Scale anchoring of one rule group, read off its rules (see Risks::assert_identifiability).
    ///
    /// The risks are kept, not counted: whether they move depends on the pinning of the run.
    /// @note Every member or pair reading the risk matches a rule of the group, asserted R-side.
    struct ScaleAnchor
    {
        bool anchored = false;                   /// Whether a rule of the group names the reference.
        std::vector<size_t> scaling_params = {}; /// Distinct risks the group applies.

        /// Whether the overall scale of the rule group is free.
        /// * a rule names the reference -> not free;
        /// * some risk applied is pinned -> not free, it fixes the level of the whole group.
        bool is_free(const Parameter &param) const;
    };

    /// Count the rule groups in `anchors` whose overall scale is free.
    size_t count_scale_free(const std::vector<ScaleAnchor> &anchors,
                            const Parameter &param);
}
