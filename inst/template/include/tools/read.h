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
 * tools/read.h
 *
 */

#ifndef HOUSEPI__TOOLS__READ__H
#define HOUSEPI__TOOLS__READ__H

#include "tools/types.h"

namespace hsp
{
    namespace read
    {
        constexpr char DELIMITER = ',';                /// CSV file delimiter.
        constexpr char LIST_DELIMITER = ':';           /// List conventional delimiter
        constexpr const char *CHAIN_COLUMN = ".chain"; /// Chain column name in draws.

        /// Extract keys from unordered_map<string, T>.
        template <typename T>
        str_vec extract_keys(const u_dict<T> &input)
        {
            str_vec entries;
            entries.reserve(input.size());
            for (const auto &pair : input)
                entries.push_back(pair.first);
            return entries;
        }

        /// Parse a list from a string.
        str_vec parse_list(const std::string &txt);

        /// Read CSV dataframe.
        u_dict<str_vec> csv_table(const std::string &pathtofile);
        /// Get number of rows in CSV table.
        size_t csv_n_rows(const u_dict<str_vec> &table);
    }
}

#endif // HOUSEPI__TOOLS__READ__H
