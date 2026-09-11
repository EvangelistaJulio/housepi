/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/types.h"
#include "function.h"
#include "group_risk.h"
#include "household.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** ContactRisk : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct DataWide;
    struct Parameter;

    struct ContactRisk : public RiskStorage
    {
    public:
        /// Build symmetric contact risks from table definitions.
        void build(DataWide &data_wide,
                   const Parameter &param,
                   const str_table &input);

        /// Whether household members i and j have matching contact risks across all members.
        /// Self-contacts are included: all contacts i<->j must also be equal.
        bool same_contacts(const Household &house,
                           MemberIdx i,
                           MemberIdx j) const noexcept
        {
            if (is_empty())
                return true;
            const Attribute &intra = at(house, i, i);
            if (intra != at(house, j, j) ||
                intra != at(house, i, j) ||
                intra != at(house, j, i))
                return false;
            for (MemberIdx k = 0; k < house.size(); ++k)
            {
                if (k == j || k == i)
                    continue;
                if (at(house, i, k) != at(house, j, k) ||
                    at(house, k, i) != at(house, k, j))
                    return false;
            }
            return true;
        }

        /// Contact risk of an ordered pair, row-major within the household.
        const Attribute &at(const Household &house,
                            MemberIdx from,
                            MemberIdx to) const noexcept
        {
            return cell_(Index::of(house),
                         static_cast<size_t>(from) * house.size() +
                             static_cast<size_t>(to));
        }

        /// Contact risk 'from' household member 'to' member.
        /// NOTE: self-contacts are set to 0.
        double value(const Household &house,
                     MemberIdx from,
                     MemberIdx to) const
        {
            return from == to ? 0. : (is_empty() ? 1. : at(house, from, to).value());
        }

        /// Contact risk intra member partition (self-contact).
        double intra_value(const Household &house, MemberIdx i) const
        {
            return is_empty() ? 1. : at(house, i, i).value();
        }

        /// Rule groups whose overall scale is free (see Parameter::assert_risk_identifiability).
        size_t n_scale_free(const Parameter &param) const
        {
            return count_scale_free_(scale_anchors_, param);
        }

    private:
        std::vector<ScaleAnchor> scale_anchors_ = {}; /// rule groups scaling the contacts
    };
}
