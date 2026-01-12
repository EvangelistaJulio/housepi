/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__GROUP__RATE__H
#define HOUSEPI__GROUP__RATE__H

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
    str_vec headers;

    GroupRate(std::shared_ptr<DataWide> dataptr,
              std::shared_ptr<Parameter> paramptr);
    ~GroupRate();

    void build(const u_dict<str_vec> &deftable);

private:
    const str_vec file_headers_ = {"group_rate", "param_name"};
    std::shared_ptr<DataWide> pDataWide_;   /// pointer to external DataWide struct
    std::shared_ptr<Parameter> pParameter_; /// pointer to external Parameter struct

    str_vec helper_assert_definitions_(const u_dict<str_vec> &deftable);
    void helper_init_all_rates_(Individual &member, const u_dict<str_vec> &these_rates);
};

#endif // HOUSEPI__GROUP__RATE__H
