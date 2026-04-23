/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "tools/read.h"
#include "tools/string.h"
#include "tools/winner.h"

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
    namespace
    {
        using idx_vec = std::vector<size_t>;

        /// Auxiliary helper to extract contact rate definition.
        struct ContactDefinitions_
        {
        public:
            const str_table &input; /// user tidy definition table
            const size_t n_rows;    /// rows in input

            std::vector<std::string> from; /// 'from'-prefixed covariates
            std::vector<std::string> to;   /// 'to'-prefixed covariates
            std::vector<std::string> data; /// covariates
            idx_vec keys;                  /// data-wide keys

            idx_vec row_gid;  /// rule group id
            idx_vec row_spec; /// rule group specificity

            /// Extract from/to definitions from input headers.
            // NOTE: 'from' and 'to' expected to have the same size
            ContactDefinitions_(
                const str_table &input_,
                const DataWide &data_wide)
                : input(input_),
                  n_rows(read::csv_n_rows(input_)),
                  row_gid(n_rows),
                  row_spec(n_rows)
            {
                const std::vector<std::string> input_headers = read::extract_keys(input);
                str::assert_no_missing_entries(input_headers, {"param_name", "rule_group"});

                for (const auto &key : input_headers)
                    if (str::starts_with(key, "from_"))
                    {
                        from.push_back(key);
                        // extract data column from 'from_' headers
                        data.push_back(key.substr(5));
                    }
                    else if (str::starts_with(key, "to_"))
                        to.push_back(key);

                helper_set_rule_groups_();
                // internally store data keys
                keys = data_wide.khold().from_string(data);
            }

            /// Associate parameter pointers to pair.
            std::vector<const Attribute *> associate_param_ptrs(
                const Parameter &param,
                const Individual &member_from,
                const Individual &member_to) const
            {
                // per-group: keep the most-specific matching rule
                std::map<size_t, Winner> best;
                for (size_t r = 0; r < n_rows; ++r)
                    if (matches_(r, member_from, member_to))
                        update_winner(best[row_gid[r]], row_spec[r], r);

                // across groups: multiply the winning params (skip reference)
                std::vector<const Attribute *> param_ptrs;
                param_ptrs.reserve(best.size());
                for (auto &[gid, w] : best)
                {
                    const std::string err_text = "contact rate from " +
                                                 member_from.err_name() + " to " +
                                                 member_to.err_name();
                    const size_t r = resolve_winner(w, input, err_text);
                    const auto &pstr = input.at("param_name")[r];
                    if (!na::is(pstr) && pstr != rates::REFERENCE)
                        param_ptrs.push_back(
                            param.at(param.khold().from_string(pstr)).pointer());
                }
                return param_ptrs;
            }

        private:
            /// Precompute per-row rule group and specificity (count of non-NA cells).
            void helper_set_rule_groups_()
            {
                for (size_t r = 0; r < n_rows; ++r)
                {
                    row_gid[r] = static_cast<size_t>(str::toi(input.at("rule_group")[r]));
                    for (size_t d = 0; d < data.size(); ++d)
                    {
                        row_spec[r] += !na::is(input.at(from[d])[r]);
                        row_spec[r] += !na::is(input.at(to[d])[r]);
                    }
                }
            }

            /// Wheter the pair of members match the contact rule at row 'r'.
            bool matches_(
                size_t r,
                const Individual &member_from,
                const Individual &member_to) const
            {
                for (size_t i = 0; i < from.size(); ++i)
                {
                    const std::string &fval = input.at(from[i])[r];
                    if (!na::is(fval) && member_from.data(keys[i]) != str::tod(fval))
                        return false;

                    const std::string &tval = input.at(to[i])[r];
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
