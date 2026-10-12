/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>
#include <vector>

#include "tools/dictionary.h"
#include "data/household.h"
#include "risks/anchor.h"
#include "risks/store.h"

//////////////////////////////  risks::Contact  //////////////////////////////

namespace housepi
{
    namespace data
    {
        struct Wide;
    }
    struct Parameter;
}

namespace housepi::risks
{
    /// Contact risks of every ordered member pair, one row-major matrix per household.
    struct Contact : public Store
    {
    public:
        //////  Setter  //////

        /// Build contact risks from rule table `input`.
        /// * contact risks are not symmetric a priori.
        void build(data::Wide &data_wide,
                   const Parameter &param,
                   const str_table &input);

        //////  Getter  //////

        /// Whether members `i` and `j` of `house` have matching contact risks with every member.
        /// * self-contacts and the contacts `i` to `j` and `j` to `i` must match too.
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

        /// Contact-risk attribute from member `from` to member `to`, row-major within `house`.
        const Attribute &at(const Household &house,
                            MemberIdx from,
                            MemberIdx to) const noexcept
        {
            return cell_(Index::of(house),
                         static_cast<size_t>(from) * house.size() +
                             static_cast<size_t>(to));
        }

        /// Contact-risk value within the partition of member `i` (self-contact attribute).
        /// * no definition -> 1.
        double intra_value(const Household &house,
                           MemberIdx i) const
        {
            return is_empty() ? 1. : at(house, i, i).value();
        }

        /// Contact-risk value from member `from` to member `to`.
        /// * self-contact -> 0;
        /// * no definition -> 1.
        double value(const Household &house,
                     MemberIdx from,
                     MemberIdx to) const
        {
            return from == to ? 0. : (is_empty() ? 1. : at(house, from, to).value());
        }

        /// Number of rule groups whose overall scale is free.
        size_t n_scale_free(const Parameter &param) const
        {
            return count_scale_free(scale_anchors_, param);
        }

    private:
        std::vector<ScaleAnchor> scale_anchors_ = {}; /// Rule groups scaling the contacts.
    };
}
