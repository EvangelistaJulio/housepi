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

        /// Parse a list from a string.
        str_vec parse_list(const std::string &txt);

        // Read CSV list.
        str_vec list_headers(const u_dict<std::string> &deflist);
        u_dict<std::string> csv_list(const std::string &pathtofile);
        // Read CSV dataframe.
        str_vec table_headers(const u_dict<str_vec> &deftable);
        u_dict<str_vec> csv_table(const std::string &pathtofile);
    }
}

#endif // HOUSEPI__TOOLS__READ__H
