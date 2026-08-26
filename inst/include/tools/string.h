/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace housepi::str
{
    constexpr size_t MAX_PRINT = 5;        /// Max elements to print
    constexpr size_t MIN_DIST = 2;         /// Similarity threshold distance
    const std::string HIDDEN_PREFIX = "."; /// Prefix for hidden objects

    /// Convert to upper case.
    std::string toup(std::string s);

    /// Convert string to boolean variable. NA is assumed as FALSE.
    bool tobool(const std::string &tkn);
    /// Convert to double using stod.
    double tod(const std::string &input, const std::string &errtxt = "");
    /// Convert to int using stoi.
    int toi(const std::string &input, const std::string &errtxt = "");
    /// Convert to unsigned long using stoul.
    unsigned long toul(const std::string &input, const std::string &errtxt = "");

    /// Format a boolean value.
    std::string format(bool condition);
    /// Format a dobule value (it deals with NA).
    std::string format(double x);
    /// Format a list of string options.
    std::string format(const std::vector<std::string> &options,
                       bool add_full_stop = true,
                       size_t max_print = MAX_PRINT);

    /// Wheter a string starts with a pattern.
    bool starts_with(std::string_view s, std::string_view pattern) noexcept;
    /// Wheter a name has been tagged as hidden.
    inline bool is_hidden(std::string_view s) noexcept { return starts_with(s, HIDDEN_PREFIX); }
    /// Make list composed of distinct entries.
    std::vector<std::string> unique(const std::vector<std::string> &v);
    /// Drop hidden strings from list.
    std::vector<std::string> drop_hidden(std::vector<std::string> v);

    /// Find the closest match, or list the options when none is close enough.
    /// The accepted distance is capped by half the length of the input.
    /// NOTE: I find 'show all available options' annoying at the moment...
    std::string find_closest(std::string_view input,
                             const std::vector<std::string> &options,
                             size_t min_dist = MIN_DIST);
}
