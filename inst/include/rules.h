/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/nan.h"
#include "tools/string.h"
#include "tools/vector.h"
#include "data/individual.h"
#include "data/wide.h"

//////////////////////////////  rules::Table  //////////////////////////////

namespace housepi::rules
{
    constexpr const char *STATUS_KEY = "status"; /// Rule column carrying the infected status.

    /// View over a rule table, with rule group, specificity and cells precomputed per row.
    ///
    /// - conflicts resolve by specificity: 'generalia specialibus non derogant';
    /// - specificity of a row: number of non-NA cells among the specificity columns;
    /// - what a row applies is read by the caller, from its own column.
    struct Table
    {
        const str_table &input;                 /// User tidy definition table.
        const std::size_t n_rows;               /// Rows in the table.
        std::vector<std::size_t> row_gid;       /// Rule-group id per row.
        std::vector<std::size_t> row_spec;      /// Specificity per row.
        std::vector<std::vector<double>> cells; /// Cells per specificity column, NA for a catch-all.

        /// View `input`, parsing `spec_columns` into the cells.
        Table(const str_table &input_,
              const std::vector<std::string> &spec_columns)
            : input(input_),
              n_rows(count_rows(input_)),
              row_gid(n_rows),
              row_spec(n_rows),
              cells(spec_columns.size(), std::vector<double>(n_rows, na::DOUBLE))
        {
            for (std::size_t r = 0; r < n_rows; ++r)
                row_gid[r] = static_cast<std::size_t>(str::toi(input.at("rule_group")[r]));

            for (std::size_t c = 0; c < spec_columns.size(); ++c)
            {
                const std::vector<std::string> &col = input.at(spec_columns[c]);
                for (std::size_t r = 0; r < n_rows; ++r)
                    if (!na::is(col[r]))
                    {
                        cells[c][r] = str::tod(col[r]);
                        ++row_spec[r];
                    }
            }
        }

        /// Whether row `r` matches `probe`, one value per specificity column.
        /// - catch-all cell: matches any value;
        /// - NA probe value: matches a catch-all only.
        bool matches(std::size_t r, const std::vector<double> &probe) const noexcept
        {
            for (std::size_t c = 0; c < cells.size(); ++c)
                if (!na::is(cells[c][r]) && probe[c] != cells[c][r])
                    return false;
            return true;
        }

        /// Pick the most specific row matching `probe`, per key.
        /// - tie at the top specificity: lowest row.
        /// @param key_of Key of row `r`.
        /// @note Ties are asserted R-side.
        template <class Key, class KeyFn>
        std::map<Key, std::size_t> resolve(const std::vector<double> &probe, KeyFn key_of) const
        {
            std::map<Key, std::size_t> best;
            for (std::size_t r = 0; r < n_rows; ++r)
                if (matches(r, probe))
                    if (auto [it, is_new] = best.try_emplace(key_of(r), r);
                        !is_new && row_spec[r] > row_spec[it->second])
                        it->second = r;
            return best;
        }
    };

    //////////////////////////////  rules::Member  //////////////////////////////

    /// Rules read against the covariates of one member.
    ///
    /// - status: carried by the member, set per probe;
    /// - every other covariate: carried by the data.
    struct Member
    {
        std::vector<std::string> data; /// Covariates carried by the data.
        std::vector<size_t> keys;      /// Data-wide keys of the covariates.
        Table rules;                   /// Rules over the covariates, then the status if keyed.

        /// Read the rules of `input` over `covariates`.
        Member(const str_table &input,
               const data::Wide &data_wide,
               const std::vector<std::string> &covariates)
            : data(vect::difference(covariates, {STATUS_KEY})),
              keys(data_wide.consult(data)),
              rules(input, vect::concatenate(data, vect::difference(covariates, data))) {}

        /// Covariates of `indiv`, then `status` if keyed.
        std::vector<double> probe(const data::Wide &data_wide,
                                  const Individual &indiv,
                                  double status) const
        {
            std::vector<double> out(rules.cells.size(), status);
            for (size_t i = 0; i < keys.size(); ++i)
                out[i] = data_wide.raw(keys[i], indiv);
            return out;
        }
    };

    //////////////////////////////  rules::CovariateClasses  //////////////////////////////

    /// Members classed by the covariates a rule reads.
    ///
    /// - members with equal covariates resolve alike;
    /// - a class resolves once for every member repeating it.
    struct CovariateClasses
    {
        std::vector<size_t> of_indiv; /// Class of every member.
        size_t count = 0;             /// Number of distinct classes.

        /// Class every member of `data_wide` on covariates `keys`.
        CovariateClasses(const data::Wide &data_wide,
                         const std::vector<size_t> &keys)
            : of_indiv(data_wide.report_indiv.total, 0)
        {
            std::map<std::vector<double>, size_t, Comparator_> class_of;
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

    private:
        /// Order covariate vectors of equal length.
        /// - missing cell: a class of its own, it matches no rule value.
        struct Comparator_
        {
            bool operator()(const std::vector<double> &a,
                            const std::vector<double> &b) const
            {
                for (size_t i = 0; i < a.size(); ++i)
                {
                    const bool na_a = na::is(a[i]);
                    const bool na_b = na::is(b[i]);
                    if (na_a != na_b)
                        return na_b;
                    if (!na_a && a[i] != b[i])
                        return a[i] < b[i];
                }
                return false;
            }
        };
    };
}
