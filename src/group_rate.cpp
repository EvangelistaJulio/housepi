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

GroupRate::GroupRate(DataWide *dataptr,
                     const Parameter *paramptr)
    : pDataWide_(dataptr),
      pParameter_(paramptr) {}
GroupRate::~GroupRate() = default;

void GroupRate::init_registry(const std::vector<std::string> &keys)
{
    pRegistry = std::make_unique<hsp::key::Registry>(keys);
}

void GroupRate::build(const hsp::str_table &input)
{
    pRegistry->clear();

    size_t n_rows = 0;
    std::vector<std::string> data_headers = {};
    if (!input.empty())
    {
        data_headers = helper_assert_definitions_(input);
        n_rows = input.begin()->second.size();
    }
    else if (pRegistry->size() == 0)
        return;

    // loop over all household members
    for (Household &house : *pDataWide_)
        for (Individual &member : house)
        {
            hsp::str_table these_rates = {};

            for (size_t r = 0; r < n_rows; ++r)
            {
                // check if all conditions are satisfied
                bool is_this_class = true;
                for (const std::string &d : data_headers)
                    if (!hsp::na::is(input.at(d)[r]))
                        is_this_class &= member.data_vec_[pDataWide_->pRegistry->from_string(d)] == hsp::str::tod(input.at(d)[r]);

                // if it's the case, add the parameter name
                if (is_this_class)
                    these_rates[input.at("group_rate")[r]].push_back(input.at("param_name")[r]);
            }

            helper_init_all_rates_(member, these_rates);
        }
}

std::vector<std::string> GroupRate::helper_assert_definitions_(const hsp::str_table &input)
{
    const std::vector<std::string> file_headers = {"group_rate", "param_name"};

    // check that tempalte headers are present
    std::vector<std::string> data_headers = hsp::read::extract_keys(input);
    hsp::str::assert_no_missing_entries(data_headers, file_headers);

    // remove template headers to keep only data
    data_headers = data_headers - file_headers;
    // check that required data columns are present
    hsp::str::assert_no_unused_entries(data_headers, pDataWide_->pRegistry->get_all());

    // assert the group rate names are allowed
    for (const std::string &groupstr : input.at("group_rate"))
    {
        if (groupstr == hsp::rates::CONTACT || groupstr == hsp::rates::MISSING)
            throw std::invalid_argument(
                "Tried to use the reserved group rate name '" + groupstr +
                "'.\nPlease, choose a different name.");
        if (hsp::na::is(groupstr))
            throw std::invalid_argument(
                "Found a missing group rate name in the definitions.");
        pRegistry->add(groupstr);
    }

    return data_headers;
}

void GroupRate::helper_init_all_rates_(Individual &member, const hsp::str_table &these_rates)
{
    // initialize group rates to reference
    member.group_rate_vec_.clear();
    member.group_rate_vec_.reserve(pRegistry->size());
    for (size_t g = 0; g < pRegistry->size(); ++g)
    {
        member.group_rate_vec_.emplace_back(hsp::rates::VALUE_REFERENCE);
        member.group_rate_vec_.back().name(pRegistry->as_string(g));
        member.group_rate_vec_.back().set_err_type("group rate");
    }

    // assign custom definitions
    for (const auto &[groupstr, vec] : these_rates)
    {
        // extract non-reference parameter pointers
        std::vector<const Attribute *> param_pointers = {};
        for (const std::string &paramstr : vec)
        {
            if (hsp::na::is(paramstr) || paramstr == hsp::rates::MISSING)
                throw std::invalid_argument(
                    "A group rate definition '" + groupstr +
                    "' was missing for " + member.err_name() + ".");
            if (paramstr != hsp::rates::REFERENCE)
                param_pointers.emplace_back(
                    pParameter_->at(pParameter_->pRegistry->from_string(paramstr)).pointer());
        }

        if (!param_pointers.empty())
        {
            size_t g = pRegistry->from_string(groupstr);
            member.group_rate_vec_[g].init_product(param_pointers);
        }
    }
}
