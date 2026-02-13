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

/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** DataLong : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::key
{
    class Registry;
}
struct DataWide;

struct DataLong
{
public:
    std::unique_ptr<hsp::key::Registry> pRegistry; /// pointer to registry

    std::vector<std::string> headers = {};
    std::vector<std::vector<std::vector<double>>> indiv_datalist = {};

    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    DataLong(const std::vector<std::string> &keys,
             DataWide *ptr);
    ~DataLong();

    //////////////////////
    //////  SETTER  //////
    //////////////////////
    
    /// Assign data (long format) after pre-processing.
    void build(const hsp::str_table &input);

private:
    DataWide *pDataWide_; // pointer to external DataWide struct

    /// Clear previous data report.
    void helper_build_clear_data_();
    /// Assert that data headers are correct.
    void helper_build_data_headers_(const hsp::str_table &input);
};
