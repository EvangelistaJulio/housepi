/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__STRING__H
#define HOUSEPI__TOOLS__STRING__H

#include <string>
#include <vector>
#include <cstdint>

#include "tools/types.h"

namespace hsp
{
    namespace str
    {
        constexpr size_t MAX_PRINT = 5;            /// Max elements to print.
        constexpr size_t MIN_DIST = 3;             /// Similarity threshold distance.
        constexpr const char *HIDDEN_PREFIX = "."; /// Prefix for hidden objects.

        /// Convert to lower case.
        std::string tolow(std::string s);
        /// Convert to upper case.
        std::string toup(std::string s);

        // converter functions
        bool tobool(const std::string &tkn);
        double tod(const std::string &input, const std::string &errtxt = "");
        int toi(const std::string &input, const std::string &errtxt = "");
        unsigned long toul(const std::string &input, const std::string &errtxt = "");
        std::int32_t toi32(const std::string &input, const std::string &errtxt);

        /// Format a dobule value (it deals with NA).
        inline std::string format(double x) { return na::is(x) ? na::STRING : std::to_string(x); }
        /// Format a boolean value.
        inline std::string format(bool condition)
        {
#ifdef _HOUSEPI_RCPP
            return (condition ? "TRUE" : "FALSE");
#else
            return (condition ? "true" : "false");
#endif
        }

        /// Format a list of string options.
        std::string format(const str_vec &options, size_t max_print = MAX_PRINT);
        /// Format a list of types castable to string.
        template <typename T>
        std::string format(const std::vector<T> &entries, size_t max_print = MAX_PRINT)
        {
            str_vec options;
            for (const T i : entries)
                options.push_back(std::to_string(i));
            return format(options, max_print);
        }
        /// Format a 2D matrix.
        std::string format(const matrix<double> &vec);

        /// Wheter a name has been tagged as private.
        bool is_private(const std::string &s);
        /// Make list composed of distinct entries.
        str_vec unique(const str_vec &v);
        /// Drop private strings from list.
        str_vec drop_private(str_vec v);

        /// Find the closest match.
        // NOTE: I find the options show alla available options annoying at the moment...
        std::string find_closest(const std::string &input,
                                 const str_vec &options,
                                 unsigned short min_dist = MIN_DIST,
                                 bool show_avail_options = false);

        /// Check for unused entries from a vector of options and throw exception.
        void assert_no_missing_entries(const str_vec &entries,
                                       const str_vec &options,
                                       const std::string &where = "column");
        void assert_no_unused_entries(const str_vec &entries,
                                      const str_vec &options,
                                      const std::string &where = "column");
    }
}

#endif // HOUSEPI__TOOLS__STRING__H
