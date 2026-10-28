/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/dictionary.h"
#include "tools/nan.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "data/individual.h"
#include "data/wide.h"
#include "function.h"
#include "parameter.h"
#include "rates/anchor.h"
#include "rates/store.h"
#include "rules.h"

//////////////////////////////  rates::read_scale_anchors  //////////////////////////////

namespace housepi::rates
{
    constexpr double VALUE_REFERENCE = 1.0; /// Rate value of the reference.

    /// Read the units each rule resolves, NA when the table gives none or it is unknown.
    inline std::vector<double> read_row_units(const str_table &input)
    {
        std::vector<double> units(count_rows(input), na::DOUBLE);
        if (const auto it = input.find("n"); it != input.end())
            for (size_t r = 0; r < units.size(); ++r)
                if (!na::is(it->second[r]))
                    units[r] = str::tod(it->second[r]);
        return units;
    }

    /// Read the scale anchoring of every rule group.
    /// * a rule naming the reference -> anchors its group, its units those of the reference;
    /// * any other rule -> a rate its group applies, its units added to the rate's.
    /// @param row_param Parameter key per rule, NA for the reference.
    /// @param row_units Units resolved by each rule, NA if unknown.
    /// @param key_of Rule-group key of row `r`.
    template <class Key, class KeyFn>
    std::map<Key, ScaleAnchor> read_scale_anchors(const std::vector<size_t> &row_param,
                                                  const std::vector<double> &row_units,
                                                  KeyFn key_of)
    {
        std::map<Key, ScaleAnchor> anchors;
        for (size_t r = 0; r < row_param.size(); ++r)
        {
            ScaleAnchor &anchor = anchors[key_of(r)];
            const size_t p = row_param[r];
            if (na::is(p))
            {
                anchor.anchored = true;
                anchor.n_reference += row_units[r];
                continue;
            }
            auto it = std::find(anchor.scaling_params.begin(), anchor.scaling_params.end(), p);
            if (it == anchor.scaling_params.end())
            {
                anchor.scaling_params.push_back(p);
                anchor.n_params.push_back(0.);
                it = anchor.scaling_params.end() - 1;
            }
            anchor.n_params[it - anchor.scaling_params.begin()] += row_units[r];
        }
        return anchors;
    }

    //////////////////////////////  rates::GroupDefinitions  //////////////////////////////

    /// Group-rate rule key: rate name and rule group.
    using GroupKey = std::pair<std::string, size_t>;

    /// Group-rate rules, read against member covariates.
    ///
    /// * status cell -> rewritten to the slot of its level.
    struct GroupDefinitions : rules::Member
    {
        std::vector<size_t> row_param = {}; /// Parameter key per rule, NA for the reference.

        /// Read the rules of `input`.
        /// @param spec_headers Every covariate a rule may key on.
        /// @param slot_of Slot of each status level.
        GroupDefinitions(const str_table &input,
                         const data::Wide &data_wide,
                         const Parameter &param,
                         const std::vector<std::string> &spec_headers,
                         const StatusTable<StatusSlot> &slot_of)
            : rules::Member(input, data_wide, spec_headers)
        {
            if (input.count("param_name"))
                row_param = param.find(input.at("param_name"));

            // a status cell reads as the slot of its level
            if (input.count(rules::STATUS_KEY))
                for (double &cell : rules.cells.back())
                    if (!na::is(cell))
                        cell = slot_of[static_cast<StatusLvl>(cell)];
        }

        /// Resolve the most specific matching rule per rate name and rule group.
        /// * covariates of `indiv` read at status slot `slot`.
        std::map<GroupKey, size_t> resolve(const data::Wide &data_wide,
                                           const Individual &indiv,
                                           StatusSlot slot) const
        {
            return rules.resolve<GroupKey>(
                probe(data_wide, indiv, slot),
                [&](size_t r)
                { return GroupKey{rules.input.at("group_rate")[r], rules.row_gid[r]}; });
        }
    };

    //////////////////////////////  rates::ContactDefinitions  //////////////////////////////

    /// Contact-rate rules, read against the covariates of an ordered member pair.
    ///
    /// * `from_` columns: covariates of the member transmitting;
    /// * `to_` columns: covariates of the member receiving, in the same order.
    struct ContactDefinitions
    {
        std::vector<std::string> data; /// Covariates, prefix dropped.
        std::vector<size_t> keys;      /// Data-wide keys of the covariates.
        rules::Table rules;            /// Rules over the `from_` columns, then the `to_` columns.
        std::vector<size_t> row_param; /// Parameter key per rule, NA for the reference.

        /// Read the rules of `input`.
        /// @note Expects as many `from_` as `to_` columns.
        ContactDefinitions(const str_table &input,
                           const data::Wide &data_wide,
                           const Parameter &param)
            : data(prefixed_(input, "from_", true)),
              keys(data_wide.consult(data)),
              rules(input, vect::concatenate(prefixed_(input, "from_"), prefixed_(input, "to_"))),
              row_param(param.find(input.at("param_name"))) {}

        /// Resolve the parameters of the rates the pair `indiv_from` to `indiv_to` multiplies.
        /// * per rule group -> the most specific matching rule;
        /// * across groups -> the product, the reference left out.
        std::vector<size_t> associate_params(const data::Wide &data_wide,
                                             const Individual &indiv_from,
                                             const Individual &indiv_to) const
        {
            std::vector<double> probe(2 * keys.size());
            for (size_t i = 0; i < keys.size(); ++i)
            {
                probe[i] = data_wide.raw(keys[i], indiv_from);
                probe[keys.size() + i] = data_wide.raw(keys[i], indiv_to);
            }
            std::vector<size_t> params;
            for (auto &[gid, r] : rules.resolve<size_t>(probe, [&](size_t row)
                                                        { return rules.row_gid[row]; }))
                if (const size_t p = row_param[r]; !na::is(p))
                    params.push_back(p);
            return params;
        }

    private:
        /// Columns of `input` starting with `prefix`, the prefix dropped when `strip`.
        static std::vector<std::string> prefixed_(const str_table &input,
                                                  std::string_view prefix,
                                                  bool strip = false)
        {
            std::vector<std::string> out;
            for (const auto &key : extract_keys(input))
                if (str::starts_with(key, prefix))
                    out.push_back(strip ? key.substr(prefix.size()) : key);
            return out;
        }
    };

    /// Rates an ordered pair of covariate classes multiplies, resolved on first encounter.
    struct PairRate
    {
        std::vector<size_t> params = {};                /// Parameter keys of the rates the pair applies.
        std::vector<const Attribute *> param_ptrs = {}; /// Rates the pair applies.
        bool resolved = false;                          /// Whether the pair was resolved.
    };
}
