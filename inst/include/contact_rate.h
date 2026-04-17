/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <vector>

#include "tools/types.h"
#include "house_storage.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** ContactRate : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct DataWide;
    struct Parameter;

    struct ContactRate
    {
    public:
        /// Build symmetric contact rates from table definitions.
        void build(DataWide &data_wide, const Parameter &param, const str_table &input);

        /// Whether no contact rate definitions were provided in the last build.
        bool is_empty() const noexcept { return house_contact_rates_.empty(); }

    private:
        std::vector<House::Contact> house_contact_rates_; /// Storage for household contact rates
    };
}
