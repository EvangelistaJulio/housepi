/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/types.h"

#include "house_storage.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** ContactRate : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
    struct DataWide;
    struct Parameter;
    class Household;
    class Individual;

    struct ContactRate
    {
    public:
        ContactRate(DataWide *dataptr, const Parameter *paramptr);
        ~ContactRate();

        /// Build symmetric contact rates from table definitions.
        void build(const str_table &input);

    private:
        std::vector<House::Contact> house_contact_rates_; /// External storage for household contact rates

        DataWide *pDataWide_;         /// pointer to external DataWide struct
        const Parameter *pParameter_; /// pointer to external Parameter struct
    };
}
