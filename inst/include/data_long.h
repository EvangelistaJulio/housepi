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

#include "indiv_storage.h"

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** DataLong : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
    namespace key
    {
        struct Registry;
    }
    struct DataWide;

    struct DataLong
    {
    public:
        std::unique_ptr<key::Registry> pRegistry; /// pointer to registry

        std::vector<std::string> headers = {};
        std::vector<Indiv::Long> indiv_datalist = {};

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        DataLong(DataWide *ptr);
        ~DataLong();

        void init_registry(const std::vector<std::string> &keys);

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign data (long format) after pre-processing.
        void build(const str_table &input);

    private:
        DataWide *pDataWide_; // pointer to external DataWide struct

        /// Clear previous data report.
        void helper_build_clear_data_();
        /// Assert that data headers are correct.
        void helper_build_data_headers_(const str_table &input);
    };
}
