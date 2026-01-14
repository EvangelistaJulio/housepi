/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__READ__H
#define HOUSEPI__TOOLS__READ__H

#include <string>
#include <vector>

#include "tools/types.h"

namespace hsp::read
{
    constexpr char DELIMITER = ',';      /// CSV file delimiter.
    constexpr char LIST_DELIMITER = ':'; /// List conventional delimiter

    /// Parse a list from a string.
    std::vector<std::string> parse_list(const std::string &txt);

    /// Get number of rows in CSV table.
    size_t csv_n_rows(const str_table &table);
    /// Read CSV dataframe.
    str_table csv_table(const std::string &pathtofile);
}

#endif // HOUSEPI__TOOLS__READ__H
