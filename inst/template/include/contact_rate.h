/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** ContactRate : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

struct DataWide;
struct Parameter;
class Household;
class Individual;

/// UPDATE: callback will be dropped in favour of table definitions
struct ContactRate
{
public:
    ContactRate(
        std::function<std::string(Household &, const Individual &, const Individual &)> func_init_contact_rate,
        DataWide *dataptr,
        const Parameter *paramptr);
    ~ContactRate();

    /// Build all group rates and contact rates and set rate per group.
    void build();

private:
    std::function<std::string(Household &, const Individual &, const Individual &)> callback_init_contact_rate_;

    DataWide *pDataWide_;         /// pointer to external DataWide struct
    const Parameter *pParameter_; /// pointer to external Parameter struct

    void init_contact_rates_(Household &house, const Individual &member_k, const Individual &member_i);
};
