/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "tools/keys.h"
#include "tools/types.h"

#include "individual.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** GroupRisk : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct DataWide;
    struct Parameter;

    /// Scale anchoring of a single rule group (see Parameter::assert_risk_identifiability).
    /// Resolving the group against a member, or an ordered pair for contact risks, either
    /// matches a rule, and reads the risk it names, or matches none, and reads the reference.
    /// The risks read are kept rather than counted: whether they can move is answered
    /// against the pinning of the run being asserted.
    struct ScaleAnchor
    {
        bool anchored = false;                   /// some resolution reads the reference
        std::vector<size_t> scaling_params = {}; /// distinct risks the group applies

        /// Register the risk read by one resolution. An empty name means no rule matched.
        void update(std::string_view pstr, const Parameter &param);
        /// Whether the overall scale of the rule group is left free:
        /// no resolution reads the reference and every risk applied can move.
        /// A single fixed risk pins the level of the whole group.
        bool is_free(const Parameter &param) const;
    };

    struct GroupRisk : public KeyModule
    {
    public:
        /// Update model report and set group risks associate to new individual.
        /// NOTE: this has to run even if input is empty to initialize existing risks to reference.
        void build(DataWide &data_wide, const Parameter &param, const str_table &deftable);
        /// Finalize all variables (see Variable::finalize).
        void finalize();

        /// Access the group risk attribute of given key for an individual.
        const Attribute &at(size_t g, const Individual &indiv) const noexcept { return unique_groups_[group_index_[Index::of(indiv)]].vec[g]; }
        /// Access the group risk value of given key for an individual.
        double value(size_t g, const Individual &indiv) const { return at(g, indiv).value(); }
        /// Whether two members belong to the same group-risk partition.
        bool same_group(const Individual &a, const Individual &b) const noexcept { return group_index_[Index::of(a)] == group_index_[Index::of(b)]; }

        /// Whether group risk storage was allocated.
        bool is_allocated() const noexcept { return !unique_groups_.empty(); }

        /// Rule groups of given key whose overall scale is free (see Parameter::assert_risk_identifiability).
        size_t n_scale_free(size_t g, const Parameter &param) const;

    private:
        std::vector<std::vector<ScaleAnchor>> scale_anchors_ = {}; /// rule groups per group risk
        /// Individual group risks.
        struct IndivGroup_
        {
            std::vector<Attribute> vec = {};
        };
        std::vector<IndivGroup_> unique_groups_; /// Unique group risk definitions
        std::vector<size_t> group_index_ = {};   /// Map individual index to unique group risk index

        /// Verify that group risks definitions are correct and pre-process theme.
        std::vector<std::string> helper_process_definitions_(DataWide &data_wide, const str_table &input);
        /// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
        void helper_init_unique_risk_(const DataWide &data_wide, const Parameter &param, const Individual &indiv, const str_table &these_risks);
    };
}
