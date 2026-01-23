/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>
#include <cstdint>

namespace hsp::str
{
    constexpr size_t MAX_PRINT = 5;                 /// Max elements to print.
    constexpr size_t MIN_DIST = 3;                  /// Similarity threshold distance.
    constexpr std::string_view HIDDEN_PREFIX = "."; /// Prefix for hidden objects.

    /// Convert to lower case.
    std::string tolow(std::string s);
    /// Convert to upper case.
    std::string toup(std::string s);

    // converter functions
    bool tobool(const std::string &tkn);
    double tod(const std::string &input, const std::string &errtxt = "");
    int toi(const std::string &input, const std::string &errtxt = "");
    unsigned long toul(const std::string &input, const std::string &errtxt = "");
    std::int32_t toi32(const std::string &input, const std::string &errtxt = "");

    /// Format a boolean value.
    std::string format(bool condition);
    /// Format a dobule value (it deals with NA).
    std::string format(double x);
    /// Format a list of string options.
    std::string format(const std::vector<std::string> &options, size_t max_print = MAX_PRINT);

    /// Wheter a name has been tagged as private.
    bool is_private(std::string_view s) noexcept;
    /// Make list composed of distinct entries.
    std::vector<std::string> unique(const std::vector<std::string> &v);
    /// Drop private strings from list.
    std::vector<std::string> drop_private(std::vector<std::string> v);

    /// Find the closest match.
    /// NOTE: I find the options show all available options annoying at the moment...
    std::string find_closest(std::string_view input,
                             const std::vector<std::string> &options,
                             unsigned short min_dist = MIN_DIST,
                             bool show_avail_options = false);

    /// Check for unused entries from a vector of options and throw exception.
    void assert_no_missing_entries(const std::vector<std::string> &entries,
                                   const std::vector<std::string> &options,
                                   const std::string &where = "column");
    void assert_no_unused_entries(const std::vector<std::string> &entries,
                                  const std::vector<std::string> &options,
                                  const std::string &where = "column");
}
