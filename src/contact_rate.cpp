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
#include "function.h"
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

void ContactRate::build()
{
    try
    {
        /// FIXME: skipped for the momemt
        return;

        // loop over all household pairs
        for (Household &house : *pDataWide_)
        {
            // initialize contact rate matrix
            const size_t size = house.size();
            house.contact_rates_.clear();
            house.contact_rates_.assign(size, std::vector<Attribute>(size, hsp::rates::VALUE_REFERENCE));
                
            for (size_t k = 0; k < size; ++k)  // from : member k
                for (size_t i = 0; i < k; ++i) // -> to : member i
                    init_contact_rates_(house, k, i);
        }
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While initializing contact rates:\n" +
                        std::string(e.what()));
    }
}

void ContactRate::init_contact_rates_(Household &house, size_t k, size_t i)
{
    const Individual &member_k = house.member(k);
    const Individual &member_i = house.member(i);
    const std::string rate = callback_init_contact_rate_(house, member_k, member_i);

    if (hsp::na::is(rate))
        return;
    else if (rate == hsp::rates::MISSING)
        throw std::invalid_argument(
            "Missing contact rate for members '" +  member_k.id() + 
            "' and '" + member_i.id() + "', from " + house.err_name() + ".");

    // default reference rate
    if (rate != hsp::rates::REFERENCE) // parameter rate
    {
        size_t p = pParameter_->pRegistry->from_string(rate);
        Attribute attr = pParameter_->at(p);
        house.contact_rates_[k][i] = attr;
        house.contact_rates_[i][k] = attr; /// NOTE: assuming symmetric
    }
}
