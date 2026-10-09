/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/dictionary.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/rule_table.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "data/individual.h"
#include "data/wide.h"
#include "function.h"
#include "parameter.h"
#include "risks/anchor.h"
#include "risks/store.h"

//////////////////////////////  risks::read_scale_anchors  //////////////////////////////

namespace housepi::risks
{
    constexpr const char *REFERENCE = ".reference"; /// Rule value choosing the reference.
    constexpr double VALUE_REFERENCE = 1.0;         /// Risk value of the reference.

    /// Read the scale anchoring of every rule group.
    /// - a rule naming the reference: anchors its group;
    /// - any other rule: a risk its group applies.
    /// @param row_param Parameter key per rule, NA for the reference.
    /// @param key_of Rule-group key of row `r`.
    template <class Key, class KeyFn>
    std::map<Key, ScaleAnchor> read_scale_anchors(const std::vector<size_t> &row_param, KeyFn key_of)
    {
        std::map<Key, ScaleAnchor> anchors;
        for (size_t r = 0; r < row_param.size(); ++r)
        {
            ScaleAnchor &anchor = anchors[key_of(r)];
            const size_t p = row_param[r];
            if (na::is(p))
                anchor.anchored = true;
            else if (!vect::contains(anchor.scaling_params, p))
                anchor.scaling_params.push_back(p);
        }
        return anchors;
    }

    //////////////////////////////  risks::GroupDefinitions  //////////////////////////////

    /// Group-risk rule key: risk name and rule group.
    using GroupKey = std::pair<std::string, size_t>;

    /// Group-risk rules, read against member covariates.
    ///
    /// - status: carried by the member, read at a status slot;
    /// - every other covariate: carried by the data.
    struct GroupDefinitions
    {
        std::vector<std::string> data;      /// Covariates carried by the data.
        std::vector<size_t> keys;           /// Data-wide keys of the covariates.
        RuleTable rules;                    /// Rules over the covariates, then the status slot if keyed.
        std::vector<size_t> row_param = {}; /// Parameter key per rule, NA for the reference.

        /// Read the rules of `input`.
        /// @param spec_headers Every covariate a rule may key on.
        /// @param slot_of Slot of each status level.
        GroupDefinitions(const str_table &input,
                         const data::Wide &data_wide,
                         const Parameter &param,
                         const std::vector<std::string> &spec_headers,
                         const StatusTable<StatusSlot> &slot_of)
            : data(vect::difference(spec_headers, {STATUS_KEY})),
              keys(data_wide.consult(data)),
              rules(input,
                    input.count(STATUS_KEY)
                        ? vect::concatenate(data, std::vector<std::string>{STATUS_KEY})
                        : data,
                    "param_name")
        {
            if (input.count("param_name"))
                row_param = param.find(input.at("param_name"));

            // a status cell reads as the slot of its level
            if (input.count(STATUS_KEY))
                for (double &cell : rules.cells.back())
                    if (!na::is(cell))
                        cell = slot_of[static_cast<StatusLvl>(cell)];
        }

        /// Resolve the most specific matching rule per risk name and rule group.
        /// - covariates of `indiv` read at status slot `slot`.
        std::map<GroupKey, size_t> resolve(const data::Wide &data_wide,
                                           const Individual &indiv,
                                           StatusSlot slot) const
        {
            std::vector<double> probe(rules.cells.size(), slot);
            for (size_t i = 0; i < keys.size(); ++i)
                probe[i] = data_wide.raw(keys[i], indiv);
            return rules.resolve<GroupKey>(
                probe,
                [&](size_t r)
                { return GroupKey{rules.input.at("group_risk")[r], rules.row_gid[r]}; },
                [&](const GroupKey &key)
                {
                    return "group risk " + msg::val(key.first) + " for " +
                           data_wide.str(indiv) + " (group " +
                           std::to_string(key.second) + ")";
                });
        }
    };

    //////////////////////////////  risks::ContactDefinitions  //////////////////////////////

    /// Contact-risk rules, read against the covariates of an ordered member pair.
    ///
    /// - `from_` columns: covariates of the member transmitting;
    /// - `to_` columns: covariates of the member receiving, in the same order.
    struct ContactDefinitions
    {
    public:
        std::vector<std::string> data; /// Covariates, prefix dropped.
        std::vector<size_t> keys;      /// Data-wide keys of the covariates.
        RuleTable rules;               /// Rules over the `from_` columns, then the `to_` columns.
        std::vector<size_t> row_param; /// Parameter key per rule, NA for the reference.

        /// Read the rules of `input`.
        /// @note Expects as many `from_` as `to_` columns.
        ContactDefinitions(const str_table &input,
                           const data::Wide &data_wide,
                           const Parameter &param)
            : data(prefixed_(input, "from_", true)),
              keys(data_wide.consult(data)),
              rules(input, vect::concatenate(prefixed_(input, "from_"), prefixed_(input, "to_")), "param_name"),
              row_param(param.find(input.at("param_name"))) {}

        /// Resolve the risks the pair `indiv_from` to `indiv_to` multiplies.
        /// - per rule group: the most specific matching rule;
        /// - across groups: the product, the reference left out.
        std::vector<const Attribute *> associate_param_ptrs(const data::Wide &data_wide,
                                                            const Parameter &param,
                                                            const Individual &indiv_from,
                                                            const Individual &indiv_to) const
        {
            std::vector<double> probe(2 * keys.size());
            for (size_t i = 0; i < keys.size(); ++i)
            {
                probe[i] = data_wide.raw(keys[i], indiv_from);
                probe[keys.size() + i] = data_wide.raw(keys[i], indiv_to);
            }
            const auto resolved = rules.resolve<size_t>(
                probe,
                [&](size_t r)
                { return rules.row_gid[r]; },
                [&](size_t)
                {
                    return "contact risk from " + data_wide.str(indiv_from) +
                           " to " + data_wide.str(indiv_to);
                });

            std::vector<const Attribute *> param_ptrs;
            for (auto &[gid, r] : resolved)
                if (const size_t p = row_param[r]; !na::is(p))
                    param_ptrs.push_back(param.at(p).pointer());
            return param_ptrs;
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

    //////////////////////////////  risks::CovariateClasses  //////////////////////////////

    /// Order covariate vectors of equal length.
    /// - missing cell: a class of its own, it matches no rule value.
    struct CovariateLess
    {
        bool operator()(const std::vector<double> &a,
                        const std::vector<double> &b) const
        {
            for (size_t i = 0; i < a.size(); ++i)
            {
                const bool na_a = na::is(a[i]), na_b = na::is(b[i]);
                if (na_a != na_b)
                    return na_b;
                if (!na_a && a[i] != b[i])
                    return a[i] < b[i];
            }
            return false;
        }
    };

    /// Members classed by the covariates a contact rule reads.
    ///
    /// - a rule reads the covariates of the two members and nothing else;
    /// - members with equal covariates resolve alike;
    /// - a pair of classes resolves once for every household repeating it.
    struct CovariateClasses
    {
    public:
        std::vector<size_t> of_indiv; /// Class of every member.
        size_t count = 0;             /// Number of distinct classes.

        /// Class every member of `data_wide` on covariates `keys`.
        CovariateClasses(const data::Wide &data_wide,
                         const std::vector<size_t> &keys)
            : of_indiv(data_wide.report_indiv.total, 0)
        {
            std::map<std::vector<double>, size_t, CovariateLess> class_of;
            std::vector<double> covariates(keys.size());
            for (const Household &house : data_wide)
                for (const Individual &indiv : house)
                {
                    for (size_t i = 0; i < keys.size(); ++i)
                        covariates[i] = data_wide.raw(keys[i], indiv);
                    of_indiv[Index::of(indiv)] =
                        class_of.emplace(covariates, class_of.size()).first->second;
                }
            count = class_of.size();
        }
    };

    /// Risks an ordered pair of covariate classes multiplies, resolved on first encounter.
    struct PairRisk
    {
        std::vector<const Attribute *> param_ptrs = {}; /// Risks the pair applies.
        bool resolved = false;                          /// Whether the pair was resolved.
    };
}
