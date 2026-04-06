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
#include "tools/operators.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

#include "group_rate.h"
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

namespace hsp
{

    GroupRate::GroupRate(DataWide *dataptr,
                         const Parameter *paramptr)
        : pDataWide_(dataptr),
          pParameter_(paramptr) {}
    GroupRate::~GroupRate() = default;

    void GroupRate::init_registry(const std::vector<std::string> &keys)
    {
        pKeys = std::make_unique<KeyHolder>(keys);
    }

    void GroupRate::build(const str_table &input)
    {
        pKeys->clear();
        indiv_group_rates_.clear();

        std::vector<std::string> data_headers = {};
        if (!input.empty())
            data_headers = helper_assert_definitions_(input);
        else if (pKeys->empty())
            return;

        // set headers to be multiplicity-relevat
        pDataWide_->set_multiplicity_data(data_headers);

        // allocate external group rate storage
        indiv_group_rates_.resize(pDataWide_->report_indiv.total);

        // loop over all household members
        size_t n_rows = read::csv_n_rows(input);
        for (Household &house : *pDataWide_)
            for (Individual &member : house)
            {
                member.pgroup_rate_ = &indiv_group_rates_[member.index_];
                str_table these_rates = {};

                for (size_t r = 0; r < n_rows; ++r)
                {
                    // check if all conditions are satisfied
                    bool match = true;
                    for (const std::string &dstr : data_headers)
                        if (!na::is(input.at(dstr)[r]))
                            match &= member.pwide_->vec[pDataWide_->pKeys->from_string(dstr)] == str::tod(input.at(dstr)[r]);

                    // if it's the case, add the parameter name
                    if (match)
                        these_rates[input.at("group_rate")[r]].push_back(input.at("param_name")[r]);
                }

                helper_init_all_rates_(member, these_rates);
            }
    }

    std::vector<std::string> GroupRate::helper_assert_definitions_(const str_table &input)
    {
        const std::vector<std::string> file_headers = {"group_rate", "param_name"};

        // check that tempalte headers are present
        std::vector<std::string> data_headers = read::extract_keys(input);
        str::assert_no_missing_entries(data_headers, file_headers);

        // remove template headers to keep only data
        data_headers = data_headers - file_headers;
        // check that required data columns are present
        str::assert_no_unused_entries(data_headers, pDataWide_->pKeys->get_all());

        // assert the group rate names are allowed
        for (const std::string &groupstr : input.at("group_rate"))
        {
            if (groupstr == rates::CONTACT || groupstr == rates::MISSING)
                throw std::invalid_argument(
                    "Tried to use the reserved group rate name '" + groupstr +
                    "'.\nPlease, choose a different name.");
            if (na::is(groupstr))
                throw std::invalid_argument(
                    "Found a missing group rate name in the definitions.");
            pKeys->add(groupstr);
        }

        return data_headers;
    }

    void GroupRate::helper_init_all_rates_(Individual &member, const str_table &these_rates)
    {
        // initialize group rates to reference (access storage directly)
        Indiv::Group &igr = indiv_group_rates_[member.index_];
        igr.vec.clear();
        igr.vec.reserve(pKeys->size());
        for (size_t g = 0; g < pKeys->size(); ++g)
        {
            igr.vec.emplace_back(rates::VALUE_REFERENCE);
            igr.vec.back().name(pKeys->as_string(g));
            igr.vec.back().set_err_type("group rate");
        }

        // assign custom definitions
        for (const auto &[groupstr, vec] : these_rates)
        {
            // extract non-reference parameter pointers
            std::vector<const Attribute *> param_pointers = {};
            for (const std::string &paramstr : vec)
            {
                if (na::is(paramstr) || paramstr == rates::MISSING)
                    throw std::invalid_argument(
                        "A group rate definition '" + groupstr +
                        "' was missing for " + member.err_name() + ".");
                if (paramstr != rates::REFERENCE)
                    param_pointers.emplace_back(
                        pParameter_->at(pParameter_->pKeys->from_string(paramstr)).pointer());
            }

            if (!param_pointers.empty())
            {
                size_t g = pKeys->from_string(groupstr);
                igr.vec[g].init_product(param_pointers);
            }
        }
    }
}
