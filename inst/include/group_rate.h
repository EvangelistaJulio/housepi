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
///////////////////////////////////  *** GroupRate : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct DataWide;
    struct Parameter;

    struct GroupRate : public KeyModule
    {
    public:
        /// Update model report and set group rates associate to new individual.
        /// NOTE: this has to run even if input is empty to initialize existing rates to reference.
        void build(DataWide &data_wide, const Parameter &param, const str_table &deftable);
        /// Optimize all variables (see Variable::optimize).
        void optimize();

        /// Access the group rate attribute of given key for a member.
        const Attribute &at(size_t g, const Individual &member) const noexcept { return unique_groups_[group_index_[Index::of(member)]].vec[g]; }
        /// Access the group rate value of given key for a member.
        double value(size_t g, const Individual &member) const { return at(g, member).value(); }
        /// Whether two members belong to the same group-rate partition.
        bool same_group(const Individual &a, const Individual &b) const noexcept { return group_index_[Index::of(a)] == group_index_[Index::of(b)]; }

    private:
        /// Individual group rates.
        struct IndivGroup_
        {
            std::vector<Attribute> vec = {};
        };
        std::vector<IndivGroup_> unique_groups_; /// Unique group rate definitions
        std::vector<size_t> group_index_ = {};   /// Map individual index to unique group rate index

        /// Verify that group rates definitions are correct and pre-process theme.
        std::vector<std::string> helper_process_definitions_(DataWide &data_wide, const str_table &input);
        /// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
        void helper_init_unique_rate_(const Parameter &param, const Individual &member, const str_table &these_rates);
    };
}
