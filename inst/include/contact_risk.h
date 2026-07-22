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

    struct ContactRisk
    {
    public:
        /// Build symmetric contact risks from table definitions.
        void build(DataWide &data_wide, const Parameter &param, const str_table &input);
        /// Finalize all variables (see Variable::finalize).
        void finalize();

        /// Whether no contact risk definitions were provided in the last build.
        bool is_empty() const noexcept { return house_contact_risks_.empty(); }

        /// Data-wide feature columns the contact-risk rules match on.
        const std::vector<std::string> &data_wide_depen() const noexcept { return data_wide_depen_; }

        /// Contact risk 'from' household member 'to' member. No self-contact (set 0 for safety).
        double value(const Household &house, size_t from, size_t to) const
        {
            if (from == to)
                return 0.;
            if (house_contact_risks_.empty())
                return 1.;
            return house_contact_risks_[Index::of(house)].mtx[from][to].value();
        }

        /// Contact risk intra member partition (self-contact).
        double intra_value(const Household &house, size_t i) const
        {
            if (house_contact_risks_.empty())
                return 1.;
            return house_contact_risks_[Index::of(house)].mtx[i][i].value();
        }

        /// Whether household members i and j have matching contact risks across all members.
        bool same_contacts(const Household &house, size_t i, size_t j) const noexcept
        {
            if (house_contact_risks_.empty())
                return true;
            const auto &mtx = house_contact_risks_[Index::of(house)].mtx;
            for (size_t k = 0; k < house.size(); ++k)
            {
                if (k == j || k == i)
                    continue;
                if (mtx[i][k] != mtx[j][k] || mtx[k][i] != mtx[k][j])
                    return false;
            }
            return true;
        }

    private:
        /// Per-household contact risk matrix.
        struct IndivContact_
        {
            std::vector<std::vector<Attribute>> mtx = {}; /// contact risks
        };
        std::vector<IndivContact_> house_contact_risks_; /// Storage for household contact risks
        std::vector<std::string> data_wide_depen_ = {};  /// data feature columns used by rules
    };
}
