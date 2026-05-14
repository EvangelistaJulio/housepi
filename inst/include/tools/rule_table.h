/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "tools/best_match.h"
#include "tools/nan.h"
#include "tools/read.h"
#include "tools/string.h"
#include "tools/types.h"

namespace housepi
{
    /// Tidy rule table: a thin view over an `str_table` carrying precomputed
    /// per-row `rule_group` id and specificity (count of non-NA feature cells).
    /// Drives specificity-based conflict resolution via best_match.h.
    struct RuleTable
    {
        const str_table &input;             /// user tidy definition table
        const std::size_t n_rows;           /// rows in input
        std::vector<std::size_t> row_gid;   /// rule_group id per row
        std::vector<std::size_t> row_spec;  /// specificity (non-NA cell count) per row

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
            std::map<Key, BestMatch> best;
            for (std::size_t r = 0; r < n_rows; ++r)
                if (matches(r))
                    update_best_match(best[key_of(r)], row_spec[r], r);

            std::map<Key, std::size_t> out;
            for (auto &[key, w] : best)
                out.emplace(key, resolve_best_match(w, input, ctx_of(key)));
            return out;
        }
    };
}
