/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/messages.h"
#include "tools/read.h"
#include "tools/string.h"

#include "contact_rate.h"
#include "function.h"
#include "data_wide.h"
#include "household.h"
#include "individual.h"
#include "parameter.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** ContactRate : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
    ContactRate::ContactRate() = default;
    ContactRate::~ContactRate() = default;

    void ContactRate::build(DataWide &data_wide, const Parameter &param, const str_table &input)
    {
        try
        {
            // UPDATE I should just allocate for restricted houses
            // reset external contact rate storage and unbind households
            house_contact_rates_.clear();
            for (size_t H = 0; H < data_wide.report_house.total; ++H)
                data_wide.at(H).pcontact_rate_ = nullptr;

            // no definitions provided -> leave empty as sentinel
            if (input.empty())
                return;

            house_contact_rates_.resize(data_wide.report_house.total);

            // parse column headers to separate from_ and to_ prefixes
            std::vector<std::string> from_headers, to_headers;
            std::vector<std::string> data_names;

            for (const auto &key : read::extract_keys(input))
            {
                if (key == "param_name")
                    continue;
                if (key.size() > 5 && key.substr(0, 5) == "from_")
                {
                    from_headers.push_back(key);
                    data_names.push_back(key.substr(5));
                }
                else if (key.size() > 3 && key.substr(0, 3) == "to_")
                {
                    to_headers.push_back(key);
                }
            }

            if (from_headers.size() != to_headers.size())
                throw std::invalid_argument(
                    "Mismatched number of 'from_' and 'to_' columns in contact rate definitions.");

            // register data columns for multiplicity
            data_wide.set_multiplicity_data(data_names);

            size_t n_rows = read::csv_n_rows(input);

            // resolve data column indices
            std::vector<size_t> data_indices = data_wide.pKeys->from_string(data_names);

            // loop over all households
            for (size_t H = 0; H < data_wide.report_house.total; ++H)
            {
                Household &house = data_wide.at(H);
                house.pcontact_rate_ = &house_contact_rates_[H];

                // initialize contact rate matrix to reference
                house_contact_rates_[H].mtx.assign(
                    house.size(), std::vector<Attribute>(house.size(), rates::VALUE_REFERENCE));

                if (n_rows == 0)
                    continue;

                for (const Individual &member_from : house)
                    for (const Individual &member_to : house)
                        if (member_from != member_to)
                            for (size_t r = 0; r < n_rows; ++r)
                            {
                                bool match = true;
                                for (size_t d = 0; d < from_headers.size(); ++d)
                                {
                                    // FIXME data indices seems here supposed equal from/to
                                    if (!na::is(input.at(from_headers[d])[r]))
                                        match &= member_from.data(data_indices[d]) ==
                                                 str::tod(input.at(from_headers[d])[r]);
                                    if (!na::is(input.at(to_headers[d])[r]))
                                        match &= member_to.data(data_indices[d]) ==
                                                 str::tod(input.at(to_headers[d])[r]);
                                }
                                if (!match)
                                    continue;

                                const std::string &paramstr = input.at("param_name")[r];
                                if (!na::is(paramstr) && paramstr != rates::REFERENCE)
                                {
                                    Attribute &attr = house_contact_rates_[H].mtx[member_from.hhindex()][member_to.hhindex()];
                                    attr.init_pointer(param.at(param.pKeys->from_string(paramstr)).pointer());
                                }
                                break;
                            }
            }
        }
        catch (const std::exception &e)
        {
            msg::error("While initializing contact rates:\n" +
                       std::string(e.what()));
        }
    }
}
