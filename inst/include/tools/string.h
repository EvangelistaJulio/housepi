/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace housepi::str
{
    constexpr size_t MAX_PRINT = 5;
    /// Similarity threshold distance.
    constexpr size_t MIN_DIST = 2;
    constexpr char LIST_DELIMITER = ':';

    /// Convert to upper case.
    std::string toup(std::string s);

    /// Convert to boolean, NA read as false.
    bool tobool(const std::string &tkn);
    /// Convert to double using `stod`.
    double tod(const std::string &input, const std::string &errtxt = "");
    /// Convert to int using `stoi`.
    int toi(const std::string &input, const std::string &errtxt = "");
    /// Convert to a 64-bit unsigned integer using `stoull`.
    std::uint64_t toull(const std::string &input, const std::string &errtxt = "");

    /// Format a boolean value.
    std::string format(bool condition);
    /// Format a double value, NA, Inf and -Inf included.
    std::string format(double x);
    /// Format a list of string options.
    std::string format(const std::vector<std::string> &options,
                       bool add_full_stop = true,
                       size_t max_print = MAX_PRINT);

    /// Whether `s` starts with `pattern`.
    bool starts_with(std::string_view s, std::string_view pattern) noexcept;
    /// Return distinct entries of `v`, in order of first appearance.
    std::vector<std::string> unique(const std::vector<std::string> &v);

    /// Suggest the closest option as a hint line.
    /// The accepted distance is capped by half the length of the input.
    /// * no option close enough -> `otherwise`, naming where the options are shown.
    std::string find_closest(std::string_view input,
                             const std::vector<std::string> &options,
                             std::string_view otherwise,
                             size_t min_dist = MIN_DIST);

    /// Remove quotes and line breaks.
    std::string clean(std::string tkn);
    /// Split a list at each `LIST_DELIMITER`, once cleaned.
    std::vector<std::string> split_list(const std::string &txt);
    /// Parse identifiers from stringified enumerators, dropping assigned values.
    std::vector<std::string> parse_names(std::string input);
}
