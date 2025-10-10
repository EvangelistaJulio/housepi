/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * tools/string.h
 *
 */

#ifndef HOUSEPI__TOOLS__STRING__H
#define HOUSEPI__TOOLS__STRING__H

#include <string>
#include <vector>

#include "tools/types.h"
#include "tools/na_utils.h"

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
#if defined(_HOUSEPI_RCPP)
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

        // make list composed of distinct entries
        str_vec unique(const str_vec &v);
        // drop private strings from list
        str_vec drop_private(str_vec v);

        // find the closest match
        std::string find_closest(const std::string &input,
                                 const str_vec &options,
                                 unsigned short min_dist = MIN_DIST,
                                 // UPDATE: I find this annoying at the moment...
                                 bool show_avail_options = false);

        // check for unused entries from a vector of options and throw exception
        void assert_no_missing_entries(const str_vec &entries,
                                       const str_vec &options,
                                       const std::string &where = "column");
        void assert_no_unused_entries(const str_vec &entries,
                                      const str_vec &options,
                                      const std::string &where = "column");
    }
}

#endif // HOUSEPI__TOOLS__STRING__H
