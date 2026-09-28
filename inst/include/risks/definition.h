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
#include "tools/string.h"
#include "tools/vector.h"
#include "data/individual.h"
#include "data/wide.h"
#include "function.h"
#include "parameter.h"
#include "risks/anchor.h"
#include "risks/rule_table.h"
#include "risks/store.h"

//////////////////////////////  risks : rule vocabulary  //////////////////////////////

namespace housepi::risks
{
    constexpr const char *REFERENCE = ".reference"; /// Rule value choosing the reference.
    constexpr const char *STATUS_KEY = "status";    /// Rule column carrying the infected status.
    constexpr double VALUE_REFERENCE = 1.0;         /// Risk value of the reference.

    /// Register on every rule group in `anchors` the risk one resolution reads.
    /// - group no rule matched: reads the reference.
    /// @param resolved Rule row each matched group resolved to.
    template <class Key>
    void update_anchors(std::map<Key, ScaleAnchor> &anchors,
                        const std::map<Key, size_t> &resolved,
                        const std::vector<std::string> &param_name,
                        const Parameter &param)
    {
        for (auto &[key, anchor] : anchors)
        {
            // empty name: no rule matched
            std::string_view pstr = {};
            if (auto found = resolved.find(key); found != resolved.end())
                pstr = param_name[found->second];
            anchor.update(pstr, param);
        }
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
    public:
        std::vector<std::string> data; /// Covariates carried by the data.
        std::vector<size_t> keys;      /// Data-wide keys of the covariates.
        RuleTable rules;               /// Rule groups and specificities.

        /// Read the rules of `input`.
        /// @param spec_headers Every covariate a rule may key on.
        /// @param slot_of Slot of each status level.
        GroupDefinitions(const str_table &input,
                         const data::Wide &data_wide,
                         const std::vector<std::string> &spec_headers,
                         const StatusTable<StatusSlot> &slot_of)
            : data(vect::difference(spec_headers, {STATUS_KEY})),
              keys(data_wide.khold(data)),
              rules(input, spec_headers),
              slot_of_(slot_of)
        {
            if (input.count(STATUS_KEY))
                status_col_ = &input.at(STATUS_KEY);
        }

        /// Resolve the most specific matching rule per risk name and rule group.
        /// - covariates of `indiv` read at status slot `slot`.
        std::map<GroupKey, size_t> resolve(const data::Wide &data_wide,
                                           const Individual &indiv,
                                           StatusSlot slot) const
        {
            return rules.resolve<GroupKey>(
                [&](size_t r)
                { return matches_(data_wide, r, indiv, slot); },
                [&](size_t r)
                { return GroupKey{rules.input.at("group_risk")[r], rules.row_gid[r]}; },
                [&](const GroupKey &key)
                {
                    return "group risk " + msg::val(key.first) + " for " +
                           data_wide.str(indiv) + " (group " +
                           std::to_string(key.second) + ")";
                });
        }

    private:
        const std::vector<std::string> *status_col_ = nullptr; /// Status column of the rules, null if absent.
        const StatusTable<StatusSlot> &slot_of_;               /// Slot each status level reads.

        /// Whether member `indiv` matches rule row `r`, read at status slot `slot`.
        bool matches_(const data::Wide &data_wide,
                      size_t r,
                      const Individual &indiv,
                      StatusSlot slot) const
        {
            if (status_col_)
            {
                const std::string &sval = (*status_col_)[r];
                if (!na::is(sval) && slot != slot_of_[str::toi(sval)])
                    return false;
            }
            for (size_t i = 0; i < data.size(); ++i)
            {
                const std::string &dval = rules.input.at(data[i])[r];
                if (!na::is(dval) && data_wide.raw(keys[i], indiv) != str::tod(dval))
                    return false;
            }
            return true;
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
        std::vector<std::string> from; /// Covariate columns prefixed `from_`.
        std::vector<std::string> to;   /// Covariate columns prefixed `to_`.
        std::vector<std::string> data; /// Covariates, prefix dropped.
        std::vector<size_t> keys;      /// Data-wide keys of the covariates.
        RuleTable rules;               /// Rule groups and specificities.

        /// Read the rules of `input`.
        /// @note Expects as many `from_` as `to_` columns.
        ContactDefinitions(const str_table &input,
                           const data::Wide &data_wide)
            : from(extract_prefixed_(input, "from_")),
              to(extract_prefixed_(input, "to_")),
              data(strip_from_prefix_(from)),
              keys(data_wide.khold(data)),
              rules(input, vect::concatenate(from, to)) {}

        /// Resolve the risks the pair `indiv_from` to `indiv_to` multiplies.
        /// - per rule group: the most specific matching rule;
        /// - across groups: the product, the reference left out.
        /// Registers the scale anchoring of every group in `anchors`.
        std::vector<const Attribute *> associate_param_ptrs(const data::Wide &data_wide,
                                                            const Parameter &param,
                                                            const Individual &indiv_from,
                                                            const Individual &indiv_to,
                                                            std::map<size_t, ScaleAnchor> &anchors) const
        {
            const auto resolved = rules.resolve<size_t>(
                [&](size_t r)
                { return matches_(data_wide, r, indiv_from, indiv_to); },
                [&](size_t r)
                { return rules.row_gid[r]; },
                [&](size_t)
                {
                    return "contact risk from " + data_wide.str(indiv_from) +
                           " to " + data_wide.str(indiv_to);
                });

            std::vector<const Attribute *> param_ptrs;
            param_ptrs.reserve(resolved.size());
            for (auto &[gid, r] : resolved)
            {
                const auto &pstr = rules.input.at("param_name")[r];
                if (!na::is(pstr) && pstr != REFERENCE)
                    param_ptrs.push_back(param.at(param.khold(pstr)).pointer());
            }

            update_anchors(anchors, resolved, rules.input.at("param_name"), param);
            return param_ptrs;
        }

    private:
        /// Columns of `input` starting with `prefix`.
        static std::vector<std::string> extract_prefixed_(const str_table &input, const char *prefix)
        {
            std::vector<std::string> out;
            for (const auto &key : extract_keys(input))
                if (str::starts_with(key, prefix))
                    out.push_back(key);
            return out;
        }
        /// Drop the `from_` prefix of every column in `from`.
        static std::vector<std::string> strip_from_prefix_(const std::vector<std::string> &from)
        {
            std::vector<std::string> out;
            out.reserve(from.size());
            for (const auto &k : from)
                out.push_back(k.substr(5));
            return out;
        }

        /// Whether the pair `indiv_from` to `indiv_to` matches rule row `r`.
        bool matches_(const data::Wide &data_wide,
                      size_t r,
                      const Individual &indiv_from,
                      const Individual &indiv_to) const
        {
            for (size_t i = 0; i < from.size(); ++i)
            {
                const std::string &fval = rules.input.at(from[i])[r];
                if (!na::is(fval) && data_wide.raw(keys[i], indiv_from) != str::tod(fval))
                    return false;

                const std::string &tval = rules.input.at(to[i])[r];
                if (!na::is(tval) && data_wide.raw(keys[i], indiv_to) != str::tod(tval))
                    return false;
            }
            return true;
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
