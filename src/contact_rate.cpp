/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/messages.h"

#include "contact_rate.h"
#include "attribute.h"
#include "data_wide.h"
#include "household.h"
#include "individual.h"
#include "parameter.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** GroupRate : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

ContactRate::ContactRate(
    std::function<std::string(Household &, const Individual &, const Individual &)> func_init_contact_rate,
    DataWide *dataptr,
    const Parameter *paramptr)
    : callback_init_contact_rate_(func_init_contact_rate),
      pDataWide_(dataptr),
      pParameter_(paramptr) {}
ContactRate::~ContactRate() = default;

/// Build all group rates and contact rates and set rate per group.
void ContactRate::build()
{
    try
    {
        // loop over all household pairs
        for (Household &house : *pDataWide_)
            for (const Individual &member_k : house)     // from : member_k
                for (const Individual &member_i : house) // -> to : member_i
                    init_contact_rates_(house, member_k, member_i);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While initializing contact rates:\n" +
                        std::string(e.what()));
    }
}

/// Set contact rate according to user assignement. // UPDATE
void ContactRate::init_contact_rates_(Household &house, const Individual &member_k, const Individual &member_i)
{
    if (member_k == member_i)
        return;

    const std::string ctc_rate = callback_init_contact_rate_(house, member_k, member_i);

    if (hsp::na::is(ctc_rate))
        return;

    // default reference rate
    Attribute y(hsp::rates::VALUE_REFERENCE);

    if (ctc_rate == hsp::rates::MISSING)
        throw std::invalid_argument(
            "Missing contact rate for individuals '" +
            member_k.id() + "' and '" + member_i.id() +
            "' (" + house.err_name() + ").");

    else if (ctc_rate != hsp::rates::REFERENCE) // parameter rate
    {
        size_t p = pParameter_->pRegistry->from_string(ctc_rate);
        y = pParameter_->at(p);
    }

    house.init_contact_rates(member_k.hhindex(), member_i.hhindex(), y);
}
