/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/operators.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"
#include "tools/best_match.h"
#include "tools/rule_table.h"

#include "contact_rate.h"
#include "group_rate.h"
#include "function.h"
#include "data_wide.h"
#include "household.h"
#include "individual.h"
#include "parameter.h"

namespace housepi
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  *** GroupRate : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    void GroupRate::build(DataWide &data_wide, const Parameter &param, const str_table &input)
    {
        // reset group rate storage and unbind individuals
        khold_.clear_optional();
        unique_groups_.clear();
        for (Household &house : data_wide)
            for (Individual &member : house)
                member.pgroup_rate_ = nullptr;

        // load required data headers or skip if no definitions provided
        std::vector<std::string> data_headers = {};
        if (!input.empty())
            data_headers = helper_process_definitions_(data_wide, input);
        else if (khold_.empty())
            return;

        // load column with group rate names, if any
        std::vector<std::string> g_input_vec = {};
        if (input.count("group_rate"))
            g_input_vec = input.at("group_rate");

        // check for undefined group rates, assumed to be equal to reference
        str_table undefined_rates = {};
        for (const std::string &gstr : khold_.get_all())
            if (std::find(g_input_vec.begin(), g_input_vec.end(), gstr) == g_input_vec.end())
                undefined_rates[gstr].push_back(rates::REFERENCE);

        // precompute rule_group / specificity and resolve column indices once
        const RuleTable rules(input, data_headers);
        const auto data_ix = data_wide.khold().from_string(data_headers);

        // predicate: does rule row r match this member's feature vector?
        using Key = std::pair<std::string, size_t>;
        auto matches_with = [&](const Individual &mb)
        {
            return [&, mb_ptr = &mb](size_t r)
            {
                for (size_t d = 0; d < data_headers.size(); ++d)
                {
                    const auto &fv = input.at(data_headers[d])[r];
                    if (!na::is(fv) && mb_ptr->pwide_->vec[data_ix[d]] != str::tod(fv))
                        return false;
                }
                return true;
            };
        };

        // unique definition map
        std::map<str_table, size_t> unique_defin_map;
        for (Household &house : data_wide)
            for (Individual &member : house)
            {
                str_table these_rates = undefined_rates;

                // per-(rate_name, rule_group): keep the most-specific matching rule
                const auto resolved = rules.resolve<Key>(
                    matches_with(member),
                    [&](size_t r) { return Key{g_input_vec[r], rules.row_gid[r]}; },
                    [&](const Key &key) {
                        return "group rate '" + key.first + "' for " +
                               member.err_name() + " (group " +
                               std::to_string(key.second) + ")";
                    });

                for (auto &[key, r] : resolved)
                    these_rates[key.first].push_back(input.at("param_name")[r]);

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
        const std::vector<std::string> file_headers = {"group_rate", "param_name", "rule_group"};

        // check that tempalte headers are present
        std::vector<std::string> data_headers = read::extract_keys(input);
        str::assert_no_missing_entries(data_headers, file_headers);
        // remove template headers to keep only data
        data_headers = data_headers - file_headers;
        // check that required data columns are present
        str::assert_no_unused_entries(data_headers, data_wide.khold().get_all());

        // assert the group rate names are allowed
        dict<size_t> definition_counts = {};
        for (const std::string &gstr : str::unique(input.at("group_rate")))
        {
            if (na::is(gstr))
                throw std::invalid_argument(
                    "Found a missing group rate name in the definitions.");
            khold_.add(gstr);
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

        // set headers to be frequency-relevant
        data_wide.set_frequency_covariates(data_headers);

        return data_headers;
    }

    void GroupRate::helper_init_unique_rate_(const Parameter &param, const Individual &member, const str_table &these_rates)
    {
        unique_groups_.emplace_back();
        Indiv::Group &igr = unique_groups_.back();
        igr.vec.reserve(khold_.size());
        for (const std::string &gstr : khold_.get_all())
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
                        param.at(param.khold().from_string(pstr)).pointer());

            if (!param_pointers.empty())
                igr.vec[khold_.from_string(gstr)].init_product(param_pointers);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** ContactRate : struct ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        using idx_vec = std::vector<size_t>;

        /// Concatenate 'from_X' and 'to_X' headers into the spec-column list.
        std::vector<std::string> concat_(
            const std::vector<std::string> &a,
            const std::vector<std::string> &b)
        {
            std::vector<std::string> out;
            out.reserve(a.size() + b.size());
            out.insert(out.end(), a.begin(), a.end());
            out.insert(out.end(), b.begin(), b.end());
            return out;
        }

        /// Auxiliary helper to extract contact rate definition.
        struct ContactDefinitions_
        {
        public:
            std::vector<std::string> from; /// 'from'-prefixed covariates
            std::vector<std::string> to;   /// 'to'-prefixed covariates
            std::vector<std::string> data; /// covariates
            idx_vec keys;                  /// data-wide keys
            RuleTable rules;               /// rule_group / specificity

            /// Extract from/to definitions from input headers.
            // NOTE: 'from' and 'to' expected to have the same size
            ContactDefinitions_(
                const str_table &input,
                const DataWide &data_wide)
                : from(extract_from_(input)),
                  to(extract_to_(input)),
                  data(strip_from_prefix_(from)),
                  keys(data_wide.khold().from_string(data)),
                  rules(validate_(input), concat_(from, to))
            {
            }

            /// Associate parameter pointers to pair.
            std::vector<const Attribute *> associate_param_ptrs(
                const Parameter &param,
                const Individual &member_from,
                const Individual &member_to) const
            {
                // per-group: keep the most-specific matching rule
                const auto resolved = rules.resolve<size_t>(
                    [&](size_t r) { return matches_(r, member_from, member_to); },
                    [&](size_t r) { return rules.row_gid[r]; },
                    [&](size_t) {
                        return "contact rate from " + member_from.err_name() +
                               " to " + member_to.err_name();
                    });

                // across groups: multiply the winning params (skip reference)
                std::vector<const Attribute *> param_ptrs;
                param_ptrs.reserve(resolved.size());
                for (auto &[gid, r] : resolved)
                {
                    const auto &pstr = rules.input.at("param_name")[r];
                    if (!na::is(pstr) && pstr != rates::REFERENCE)
                        param_ptrs.push_back(
                            param.at(param.khold().from_string(pstr)).pointer());
                }
                return param_ptrs;
            }

        private:
            static const str_table &validate_(const str_table &input)
            {
                str::assert_no_missing_entries(
                    read::extract_keys(input), {"param_name", "rule_group"});
                return input;
            }
            static std::vector<std::string> extract_from_(const str_table &input)
            {
                std::vector<std::string> out;
                for (const auto &key : read::extract_keys(input))
                    if (str::starts_with(key, "from_"))
                        out.push_back(key);
                return out;
            }
            static std::vector<std::string> extract_to_(const str_table &input)
            {
                std::vector<std::string> out;
                for (const auto &key : read::extract_keys(input))
                    if (str::starts_with(key, "to_"))
                        out.push_back(key);
                return out;
            }
            static std::vector<std::string> strip_from_prefix_(
                const std::vector<std::string> &from)
            {
                std::vector<std::string> out;
                out.reserve(from.size());
                for (const auto &k : from)
                    out.push_back(k.substr(5));
                return out;
            }

            /// Whether the pair of members match the contact rule at row 'r'.
            bool matches_(
                size_t r,
                const Individual &member_from,
                const Individual &member_to) const
            {
                for (size_t i = 0; i < from.size(); ++i)
                {
                    const std::string &fval = rules.input.at(from[i])[r];
                    if (!na::is(fval) && member_from.data(keys[i]) != str::tod(fval))
                        return false;

                    const std::string &tval = rules.input.at(to[i])[r];
                    if (!na::is(tval) && member_to.data(keys[i]) != str::tod(tval))
                        return false;
                }
                return true;
            };
        };
    }

    void ContactRate::build(
        DataWide &data_wide,
        const Parameter &param,
        const str_table &input)
    {
        // reset storage and unbind households
        house_contact_rates_.clear();
        for (size_t H = 0; H < data_wide.report_house.total; ++H)
            data_wide.at(H).pcontact_rate_ = nullptr;
        if (input.empty())
            return;

        // load contact rate definitions
        ContactDefinitions_ defs(input, data_wide);

        // set headers to be frequency-relevant
        data_wide.set_frequency_covariates(defs.data);

        // loop over all households
        // UPDATE: this should be just over distinct_head households.
        house_contact_rates_.resize(data_wide.report_house.total);
        for (size_t H = 0; H < data_wide.report_house.total; ++H)
        {
            Household &house = data_wide.at(H);
            size_t house_size = house.size();

            house.pcontact_rate_ = &house_contact_rates_[H];
            house_contact_rates_[H].mtx.assign(
                house_size, std::vector<Attribute>(house_size, rates::VALUE_REFERENCE));

            // NOTE: a priori contact rates are non-symmetric
            for (size_t from = 0; from < house_size; ++from)
                for (size_t to = 0; to < house_size; ++to)
                {
                    std::vector<const Attribute *> param_ptrs =
                        defs.associate_param_ptrs(param, house.member(from), house.member(to));

                    if (!param_ptrs.empty())
                        house_contact_rates_[H].mtx[from][to].init_product(param_ptrs);
                }
        }
    }
}
