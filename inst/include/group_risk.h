/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
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

    struct GroupRisk : public KeyModule
    {
    public:
        /// Update model report and set group risks associate to new individual.
        /// NOTE: this has to run even if input is empty to initialize existing risks to reference.
        void build(DataWide &data_wide, const Parameter &param, const str_table &deftable);
        /// Finalize all variables (see Variable::finalize).
        void finalize();

        /// Access the group risk attribute of given key for a member.
        const Attribute &at(size_t g, const Individual &member) const noexcept { return unique_groups_[group_index_[Index::of(member)]].vec[g]; }
        /// Access the group risk value of given key for a member.
        double value(size_t g, const Individual &member) const { return at(g, member).value(); }
        /// Whether two members belong to the same group-risk partition.
        bool same_group(const Individual &a, const Individual &b) const noexcept { return group_index_[Index::of(a)] == group_index_[Index::of(b)]; }

        /// Data-wide feature columns the group-risk rules match on.
        const std::vector<std::string> &data_wide_depen() const noexcept { return data_wide_depen_; }

        /// Rule groups of given key whose overall scale is free (see Parameter::assert_risk_identifiability).
        size_t n_scale_free(size_t g) const noexcept { return n_scale_free_[g]; }

    private:
        std::vector<std::string> data_wide_depen_ = {}; /// data feature columns used by rules
        std::vector<size_t> n_scale_free_ = {};               /// free-scale rule groups per group risk
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
        void helper_init_unique_risk_(const DataWide &data_wide, const Parameter &param, const Individual &member, const str_table &these_risks);
    };
}
