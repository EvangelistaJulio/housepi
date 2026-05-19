/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
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
    using BestMatch = std::pair<std::size_t, std::vector<std::size_t>>;

    /// Register a matching row under its key's current best match.
    inline void update_best_match(
        BestMatch &b,
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

    /// Return the surviving row id; throw on tie at max specificity.
    inline std::size_t resolve_best_match(
        const BestMatch &b,
        const str_table &input,
        const std::string &ctx)
    {
        if (b.second.size() == 1)
            return b.second[0];
        std::string names;
        for (std::size_t r : b.second)
            names += (names.empty() ? "" : ", ") + input.at("param_name")[r];
        throw std::invalid_argument(
            "Ambiguous " + ctx + ": tied rules [" + names + "].\n" +
            "Add a disambiguating rule in order to proceed.");
    }
}
