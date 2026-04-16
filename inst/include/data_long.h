/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/types.h"

#include "indiv_storage.h"

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** DataLong : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct DataWide;

    struct DataLong : public KeyModule
    {
    public:
        std::vector<std::string> headers = {};
        std::vector<Indiv::Long> indiv_datalist = {};

        /// Assign data (long format) after pre-processing.
        void build(DataWide &data_wide, const str_table &input);

    private:
        /// Clear previous data report.
        void helper_build_clear_data_();
        /// Assert that data headers are correct.
        void helper_build_data_headers_(const DataWide &data_wide, const str_table &input);
    };
}
