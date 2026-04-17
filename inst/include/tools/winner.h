/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "tools/types.h"

namespace housepi
{
    /// Specificity-based conflict resolution
    /// 'generalia specialibus non derogant'

    /// Per-key best rule: (max specificity, rows tied at it).
    using Winner = std::pair<std::size_t, std::vector<std::size_t>>;

    /// Register a matching row under its key's current winner.
    inline void update_winner(
        Winner &w,
        std::size_t spec,
        std::size_t row)
    {
        if (w.second.empty() || spec > w.first)
        {
            w.first = spec;
            w.second.assign(1, row);
        }
        else if (spec == w.first)
            w.second.push_back(row);
    }

    /// Return the surviving row id; throw on tie at max specificity.
    inline std::size_t resolve_winner(
        const Winner &w,
        const str_table &input,
        const std::string &ctx)
    {
        if (w.second.size() == 1)
            return w.second[0];
        std::string names;
        for (std::size_t r : w.second)
            names += (names.empty() ? "" : ", ") + input.at("param_name")[r];
        throw std::invalid_argument(
            "Ambiguous " + ctx + ": tied rules [" + names + "].\n" +
            "Add a disambiguating rule in order to proceed.");
    }
}
