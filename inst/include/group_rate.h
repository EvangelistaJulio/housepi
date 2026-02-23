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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** GroupRate : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::key
{
    class Registry;
}

struct DataWide;
struct Parameter;
class Individual;

struct GroupRate
{
public:
    std::unique_ptr<hsp::key::Registry> pRegistry; /// pointer to registry

    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////
    
    GroupRate(const std::vector<std::string> &keys,
              DataWide *dataptr,
              const Parameter *paramptr);
    ~GroupRate();

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// Update model report and set group rates associate to new individual.
    /// NOTE: this has to run even if input is empty to initialize existing rates to reference.
    void build(const hsp::str_table &deftable);

private:
    DataWide *pDataWide_;   /// pointer to external DataWide struct
    const Parameter *pParameter_; /// pointer to external Parameter struct

    /// Verify that group rates definitions are correct.
    std::vector<std::string> helper_assert_definitions_(const hsp::str_table &input);
    /// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
    void helper_init_all_rates_(Individual &member, const hsp::str_table &these_rates);
};
