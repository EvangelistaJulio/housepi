/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>

#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Translator : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    namespace parse
    {
        struct Entry;
    }
    struct DataWide;
    struct GroupRisk;
    struct Augmented;
    struct Parameter;
    class Individual;
    class Household;

    /// Cross-module translator of parsed definition into Attribute.
    struct Translator
    {
        Translator(Parameter &param,
                   Augmented &augm)
            : param_(param), augm_(augm) {}

        /// Bind data-wide storage and group-risk storage for member-specific source resolution.
        void bind_data_wide(const DataWide *data_wide) noexcept { pdata_wide_ = data_wide; }
        void bind_group_risk(const GroupRisk *group_risk) noexcept { pgroup_risk_ = group_risk; }

        /// Translate a definition node or a string to a generic attribute.
        /// Within a household 'house_ptr' and no member bound, a per-member data cell
        /// is only read through a reduction folding it over the members.
        const Attribute &to_attribute(const parse::Entry &node,
                                      const Individual *indiv_ptr = nullptr,
                                      const Household *house_ptr = nullptr) const;
        const Attribute &to_attribute(const std::string &expr,
                                      const Individual *indiv_ptr = nullptr) const;

        /// Translate a definition node or a string to a numeric constant attribute.
        const Attribute &to_numeric(const parse::Entry &node,
                                    const Individual *indiv_ptr = nullptr) const;
        const Attribute &to_numeric(const std::string &expr,
                                    const Individual *indiv_ptr = nullptr) const;

        /// If node is parameter, mark as producer of transmission cache.
        void mark_param_produces_transm_cache(const parse::Entry &node) const;

        /// Whether reduction 'node' folds over the members of household 'house_ptr'.
        bool is_house_fold(const parse::Entry &node,
                           const Household *house_ptr) const;
        /// Bind the arguments of function 'node' into 'parent'. A household fold binds
        /// each operand reading a per-member cell once per member, skipping empty ones.
        /// Return whether the call is defined within a household: every operand outside
        /// a fold carries a value, and a fold holds at least one.
        bool bind_arguments(Attribute &parent,
                            const parse::Entry &node,
                            const Individual *indiv_ptr,
                            const Household *house_ptr) const;

    private:
        Parameter &param_;             /// parameter key-module
        Augmented &augm_;              /// latent key-module
        const DataWide *pdata_wide_;   /// data-wide storage and key-holder
        const GroupRisk *pgroup_risk_; /// group-risks storage and key-holder
        mutable Attribute storage_;    /// auxiliary storage attribute

        void helper_store_read_function_(const parse::Entry &node,
                                         const Individual *indiv_ptr,
                                         const Household *house_ptr) const;
        /// Whether node reads a per-member data cell outside an inner reduction.
        bool reads_member_cell_(const parse::Entry &node) const;
        bool reads_member_cell_(const std::vector<parse::Entry> &nodes) const;
    };
}
