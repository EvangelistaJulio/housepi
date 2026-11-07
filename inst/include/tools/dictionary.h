/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace housepi
{
    /// Map keyed by name, with lookup by string view.
    template <typename T>
    using dict = std::map<std::string, T, std::less<>>;
    /// String columns, by name.
    using str_table = dict<std::vector<std::string>>;
    /// Double columns, by name.
    using dbl_table = dict<std::vector<double>>;

    /// Count rows of `table`, read from its first column.
    template <typename T>
    size_t count_rows(const dict<std::vector<T>> &table) noexcept
    {
        return table.empty() ? 0 : table.begin()->second.size();
    }

    /// Return keys of `input`, in order.
    template <typename T>
    std::vector<std::string> extract_keys(const dict<T> &input)
    {
        std::vector<std::string> keys;
        keys.reserve(input.size());
        for (const auto &[key, val] : input)
            keys.push_back(key);
        return keys;
    }
}
