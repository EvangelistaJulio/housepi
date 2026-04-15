/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/types.h"

namespace housepi::read
{
    constexpr char DELIMITER = ',';      /// CSV file delimiter.
    constexpr char LIST_DELIMITER = ':'; /// List conventional delimiter

    /// Clean string from forbidden characters and convert to UNDEFINED.
    std::string strcleaner(std::string tkn);

    /// Parse a list from a string.
    std::vector<std::string> parse_list(const std::string &txt);

    /// Parse comma-separated enum names from stringified variadic arguments.
    std::vector<std::string> helper_parse_names(std::string input);

    /// Get number of rows in CSV table.
    template <typename T>
    size_t csv_n_rows(const dict<std::vector<T>> &table) noexcept
    {
        if (table.empty())
            return 0;
        return table.begin()->second.size();
    }

    /// Extract keys from dictionary.
    template <typename T>
    std::vector<std::string> extract_keys(const dict<T> &input)
    {
        std::vector<std::string> entries;
        entries.reserve(input.size());
        for (const auto &[key, val] : input)
            entries.push_back(key);
        return entries;
    }
}
