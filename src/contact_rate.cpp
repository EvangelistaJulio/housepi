/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
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

namespace housepi
{
    ContactRate::ContactRate() = default;
    ContactRate::~ContactRate() = default;

    void ContactRate::build(DataWide &data_wide, const Parameter &param, const str_table &input)
    {
        // reset contact rate storage and unbind households
        house_contact_rates_.clear();
        for (size_t H = 0; H < data_wide.report_house.total; ++H)
            data_wide.at(H).pcontact_rate_ = nullptr;

        // skip if no definitions provided
        if (input.empty())
            return;

        // parse column headers to separate 'from_' and 'to_' prefixes
        std::vector<std::string> from_headers, to_headers;
        for (const auto &key : read::extract_keys(input))
            if (str::starts_with(key, "from_"))
                from_headers.push_back(key);
            else if (str::starts_with(key, "to_"))
                to_headers.push_back(key);
            else
                continue;
        if (from_headers.size() != to_headers.size())
            throw std::invalid_argument(
                "Mismatching number of 'from_' and 'to_' columns in contact rate definitions.");

        // extract corresponding data names
        std::vector<std::string> data_names;
        for (const std::string &from : from_headers)
            data_names.push_back(from.substr(5));

        // register data columns for multiplicity
        data_wide.set_multiplicity_data(data_names);

        // resolve data column indices
        std::vector<size_t> data_indices = data_wide.khold().from_string(data_names);

        // loop over all households
        // UPDATE: this should be just over restricted households.
        house_contact_rates_.resize(data_wide.report_house.total);
        const size_t n_rows = read::csv_n_rows(input);
        for (size_t H = 0; H < data_wide.report_house.total; ++H)
        {
            Household &house = data_wide.at(H);
            house.pcontact_rate_ = &house_contact_rates_[H];

            // initialize contact rate matrix to reference
            house_contact_rates_[H].mtx.assign(
                house.size(), std::vector<Attribute>(house.size(), rates::VALUE_REFERENCE));

            // NOTE: a priori contact rates are non-symmetric
            for (size_t from = 0; from < house.size(); ++from)
                for (size_t to = 0; to < house.size(); ++to)
                {
                    if (from == to)
                        continue;

                    const Individual &member_from = house.member(from);
                    const Individual &member_to = house.member(to);

                    std::vector<const Attribute *> param_ptrs = {};
                    for (size_t r = 0; r < n_rows; ++r)
                    {
                        bool match = true;
                        for (size_t d = 0; d < from_headers.size(); ++d)
                        {
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
                            param_ptrs.push_back(
                                param.at(param.khold().from_string(paramstr)).pointer());
                    }

                    if (!param_ptrs.empty())
                        house_contact_rates_[H].mtx[from][to].init_product(param_ptrs);
                }
        }
    }
}
