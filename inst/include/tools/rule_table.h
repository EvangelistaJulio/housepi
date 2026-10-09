/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tools/dictionary.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/string.h"

//////////////////////////////  RuleTable  //////////////////////////////

namespace housepi
{
    constexpr const char *STATUS_KEY = "status"; /// Rule column carrying the infected status.

    /// View over a rule table, with rule group, specificity and cells precomputed per row.
    ///
    /// - conflicts resolve by specificity: 'generalia specialibus non derogant';
    /// - specificity of a row: number of non-NA cells among the specificity columns;
    /// - what a row applies is read by the caller, from its own column.
    struct RuleTable
    {
        const str_table &input;                 /// User tidy definition table.
        const std::size_t n_rows;               /// Rows in the table.
        std::vector<std::size_t> row_gid;       /// Rule-group id per row.
        std::vector<std::size_t> row_spec;      /// Specificity per row.
        std::vector<std::vector<double>> cells; /// Cells per specificity column, NA for a catch-all.

        /// View `input`, parsing `spec_columns` into the cells.
        /// @param label_col Column naming a row in the tie error.
        RuleTable(const str_table &input_,
                  const std::vector<std::string> &spec_columns,
                  std::string label_col)
            : input(input_),
              n_rows(count_rows(input_)),
              row_gid(n_rows),
              row_spec(n_rows),
              cells(spec_columns.size(), std::vector<double>(n_rows, na::DOUBLE)),
              label_col_(std::move(label_col))
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
        /// - tie at the top specificity: error.
        /// @param key_of Key of row `r`.
        /// @param ctx_of Error context of key `k`, read on a tie.
        template <class Key, class KeyFn, class CtxFn>
        std::map<Key, std::size_t> resolve(const std::vector<double> &probe,
                                           KeyFn key_of,
                                           CtxFn ctx_of) const
        {
            std::map<Key, BestMatch_> best;
            for (std::size_t r = 0; r < n_rows; ++r)
                if (matches(r, probe))
                    update_best_match_(best[key_of(r)], row_spec[r], r);

            std::map<Key, std::size_t> out;
            for (auto &[key, w] : best)
                out.emplace(key, resolve_best_match_(w, ctx_of(key)));
            return out;
        }

    private:
        std::string label_col_; /// Column naming a row in the tie error.

        /// Best rule of one key: top specificity, and the rows tied at it.
        using BestMatch_ = std::pair<std::size_t, std::vector<std::size_t>>;

        /// Register matching row `row` against the current best match of its key.
        static void update_best_match_(BestMatch_ &b,
                                       std::size_t spec,
                                       std::size_t row)
        {
            if (b.second.empty() || spec > b.first)
            {
                b.first = spec;
                b.second.assign(1, row);
            }
            else if (spec == b.first)
                b.second.push_back(row);
        }

        /// Return the surviving row, throw on a tie at the top specificity.
        std::size_t resolve_best_match_(const BestMatch_ &b,
                                        const std::string &ctx) const
        {
            if (b.second.size() == 1)
                return b.second[0];
            std::string names;
            for (std::size_t r : b.second)
                names += (names.empty() ? "" : ", ") + input.at(label_col_)[r];
            msg::stop(
                "Ambiguous " + ctx + ": tied rules [" + names + "].\n" +
                "Add a disambiguating rule in order to proceed.");
        }
    };
}
