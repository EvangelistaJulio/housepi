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

#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

namespace housepi
{
    /// Thin view over `str_table`.
    /// It carries precomputed per-row `rule_group` id and specificity.
    /// Drives specificity-based conflict resolution:
    /// 'generalia specialibus non derogant'
    struct RuleTable
    {
        const str_table &input;            /// user tidy definition table
        const std::size_t n_rows;          /// rows in input
        std::vector<std::size_t> row_gid;  /// rule_group id per row
        std::vector<std::size_t> row_spec; /// specificity (non-NA cell count) per row

        /// `spec_columns`: input columns counted into the specificity score.
        RuleTable(const str_table &input_,
                  const std::vector<std::string> &spec_columns)
            : input(input_),
              n_rows(read::csv_n_rows(input_)),
              row_gid(n_rows),
              row_spec(n_rows)
        {
            for (std::size_t r = 0; r < n_rows; ++r)
            {
                row_gid[r] = static_cast<std::size_t>(
                    str::toi(input.at("rule_group")[r]));
                for (const std::string &col : spec_columns)
                    row_spec[r] += !na::is(input.at(col)[r]);
            }
        }

        /// Pick the most-specific matching row per Key.
        /// `matches(r)`: row predicate.
        /// `key_of(r)`:  Key extracted from row r.
        /// `ctx_of(k)`:  error-context builder, used on tie at max specificity.
        template <class Key, class Predicate, class KeyFn, class CtxFn>
        std::map<Key, std::size_t>
        resolve(Predicate matches, KeyFn key_of, CtxFn ctx_of) const
        {
            std::map<Key, BestMatch_> best;
            for (std::size_t r = 0; r < n_rows; ++r)
                if (matches(r))
                    update_best_match_(best[key_of(r)], row_spec[r], r);

            std::map<Key, std::size_t> out;
            for (auto &[key, w] : best)
                out.emplace(key, resolve_best_match_(w, input, ctx_of(key)));
            return out;
        }

    private:
        /// Per-key best rule: max specificity, rows tied at it.
        using BestMatch_ = std::pair<std::size_t, std::vector<std::size_t>>;

        /// Register a matching row under its key's current best match.
        static void update_best_match_(
            BestMatch_ &b,
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

        /// Return the surviving row id. Throw on tie at max specificity.
        static std::size_t resolve_best_match_(
            const BestMatch_ &b,
            const str_table &input,
            const std::string &ctx)
        {
            if (b.second.size() == 1)
                return b.second[0];
            std::string names;
            for (std::size_t r : b.second)
                names += (names.empty() ? "" : ", ") + input.at("param_name")[r];
            msg::stop(
                "Ambiguous " + ctx + ": tied rules [" + names + "].\n" +
                "Add a disambiguating rule in order to proceed.");
        }
    };
}
