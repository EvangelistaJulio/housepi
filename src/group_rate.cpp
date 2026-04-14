/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/messages.h" // FLAG
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

    GroupRate::GroupRate() = default;
    GroupRate::~GroupRate() = default;

    void GroupRate::init_registry(const std::vector<std::string> &keys)
    {
        pKeys = std::make_unique<KeyHolder>(keys);
    }

    void GroupRate::build(DataWide &data_wide, const Parameter &param, const str_table &input)
    {
        // reset group rate storage and unbind individuals
        pKeys->clear_optional();
        unique_groups_.clear();
        for (Household &house : data_wide)
            for (Individual &member : house)
                member.pgroup_rate_ = nullptr;

        // load required data headers or skip if no definitions provided
        std::vector<std::string> data_headers = {};
        if (!input.empty())
            data_headers = helper_process_definitions_(data_wide, input);
        else if (pKeys->empty())
            return;

        // load column with group rate names, if any
        std::vector<std::string> g_input_vec = {};
        if (input.count("group_rate"))
            g_input_vec = input.at("group_rate");

        // check for undefined group rates, assumed to be equal to reference
        str_table undefined_rates = {};
        for (const std::string &gstr : pKeys->get_all())
            if (std::find(g_input_vec.begin(), g_input_vec.end(), gstr) == g_input_vec.end())
                undefined_rates[gstr].push_back(rates::REFERENCE);

        // unique definition map
        std::map<str_table, size_t> unique_defin_map;
        const size_t n_rows = read::csv_n_rows(input);
        for (Household &house : data_wide)
            for (Individual &member : house)
            {
                str_table these_rates = undefined_rates;

                for (size_t r = 0; r < n_rows; ++r)
                {
                    bool match = true;
                    for (const std::string &dstr : data_headers)
                        if (!na::is(input.at(dstr)[r]))
                            match &= member.pwide_->vec[data_wide.pKeys->from_string(dstr)] == str::tod(input.at(dstr)[r]);

                    if (match)
                        these_rates[g_input_vec[r]].push_back(input.at("param_name")[r]);
                }

                auto it = unique_defin_map.find(these_rates);
                if (it == unique_defin_map.end())
                {
                    it = unique_defin_map.emplace(these_rates, unique_groups_.size()).first;
                    helper_init_unique_rate_(param, member, these_rates);
                }
                member.pgroup_rate_ = &unique_groups_[it->second];
            }
    }

    std::vector<std::string> GroupRate::helper_process_definitions_(DataWide &data_wide, const str_table &input)
    {
        const std::vector<std::string> file_headers = {"group_rate", "param_name"};

        // check that tempalte headers are present
        std::vector<std::string> data_headers = read::extract_keys(input);
        str::assert_no_missing_entries(data_headers, file_headers);
        // remove template headers to keep only data
        data_headers = data_headers - file_headers;
        // check that required data columns are present
        str::assert_no_unused_entries(data_headers, data_wide.pKeys->get_all());

        // assert the group rate names are allowed
        dict<size_t> definition_counts = {};
        for (const std::string &gstr : str::unique(input.at("group_rate")))
        {
            if (gstr == rates::CONTACT)
                throw std::invalid_argument(
                    "Tried to use the reserved group rate name '" + gstr +
                    "'.\nPlease, choose a different name.");
            if (na::is(gstr))
                throw std::invalid_argument(
                    "Found a missing group rate name in the definitions.");
            pKeys->add(gstr);
            definition_counts.emplace(gstr, 0);
        }

        // allocate unique group rate definitions
        // NOTE: reserve cartesian product of distinct definitions per group rate
        for (const std::string &gstr : input.at("group_rate"))
            ++definition_counts[gstr];
        size_t max_unique = 1;
        for (auto &[gstr, count] : definition_counts)
            max_unique *= count;
        unique_groups_.reserve(max_unique);

        // set headers to be multiplicity-relevat
        data_wide.set_multiplicity_data(data_headers);

        return data_headers;
    }

    void GroupRate::helper_init_unique_rate_(const Parameter &param, const Individual &member, const str_table &these_rates)
    {
        unique_groups_.emplace_back();
        Indiv::Group &igr = unique_groups_.back();
        igr.vec.reserve(pKeys->size());
        for (const std::string &gstr : pKeys->get_all())
        {
            if (!these_rates.count(gstr))
                throw std::invalid_argument(
                    "A group rate definition '" + gstr +
                    "' was missing for " + member.err_name() + ".");

            igr.vec.emplace_back(rates::VALUE_REFERENCE);
            igr.vec.back().name(gstr);
            igr.vec.back().set_err_type("group rate");
        }

        // assign custom definitions
        for (const auto &[gstr, pstr_vec] : these_rates)
        {
            // extract non-reference parameter pointers
            std::vector<const Attribute *> param_pointers = {};
            for (const std::string &pstr : pstr_vec)
                if (pstr != rates::REFERENCE)
                    param_pointers.emplace_back(
                        param.at(param.pKeys->from_string(pstr)).pointer());

            if (!param_pointers.empty())
                igr.vec[pKeys->from_string(gstr)].init_product(param_pointers);
        }
    }
}
