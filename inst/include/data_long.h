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
        /// User defined data-long headers
        std::vector<std::string> headers = {};

        /// Individual long-format data.
        struct IndivLong
        {
            std::vector<std::vector<double>> vec = {};
        };
        std::vector<IndivLong> indiv_datalist = {};

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign data (long format) after pre-processing.
        void build(DataWide &data_wide, const str_table &input);

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Access all long-format data for a member.
        const std::vector<std::vector<double>> *get_data_long_ptr(const Individual &member) const noexcept
        {
            return &indiv_datalist[Index::of(member)].vec;
        }

    private:
        /// Clear previous data report.
        void helper_build_clear_data_();
        /// Assert that data headers are correct.
        void helper_build_data_headers_(const DataWide &data_wide, const str_table &input);
    };
}
