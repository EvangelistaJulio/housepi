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

struct DataWide;
struct Parameter;
class Individual;

struct GroupRate
{
public:
    std::shared_ptr<hsp::key::Registry> pRegistry; /// pointer to registry

    GroupRate(const std::vector<std::string> &keys,
              std::shared_ptr<DataWide> dataptr,
              std::shared_ptr<Parameter> paramptr);
    ~GroupRate();

    void build(const hsp::str_table &deftable);

private:
    std::shared_ptr<DataWide> pDataWide_;   /// pointer to external DataWide struct
    std::shared_ptr<Parameter> pParameter_; /// pointer to external Parameter struct

    std::vector<std::string> helper_assert_definitions_(const hsp::str_table &input);
    void helper_init_all_rates_(Individual &member, const hsp::str_table &these_rates);
};
