/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <string_view>
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

    /// Scale anchoring of one rule group (see Risks::assert_identifiability).
    ///
    /// Resolving the group against a member, or an ordered pair for contact risks:
    /// - matches a rule: reads the risk it names;
    /// - matches none: reads the reference.
    /// The risks read are kept, not counted: whether they move depends on the pinning of the run.
    struct ScaleAnchor
    {
        bool anchored = false;                   /// Whether some resolution reads the reference.
        std::vector<size_t> scaling_params = {}; /// Distinct risks the group applies.

        /// Register the risk one resolution reads.
        /// @param pstr Parameter name of the matched rule, empty when no rule matched.
        void update(std::string_view pstr,
                    const Parameter &param);
        /// Whether the overall scale of the rule group is free.
        /// - some resolution reads the reference: not free;
        /// - some risk applied is pinned: not free, it fixes the level of the whole group.
        bool is_free(const Parameter &param) const;
    };

    /// Count the rule groups in `anchors` whose overall scale is free.
    size_t count_scale_free(const std::vector<ScaleAnchor> &anchors,
                            const Parameter &param);
}
