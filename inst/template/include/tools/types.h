/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__TYPES__H
#define HOUSEPI__TOOLS__TYPES__H

/// DEVELOPER: flag used to develop Rcpp chunks
//#define _HOUSEPI_RCPP

#include <string>
#include <vector>
#include <unordered_map>

namespace hsp
{
    using named_args = std::vector<std::pair<std::string, double>>;
    template <typename T>
    using dict = std::unordered_map<std::string, T>;
    using str_table = dict<std::vector<std::string>>;
    using dbl_table = dict<std::vector<double>>;
}

#endif // HOUSEPI__TOOLS__TYPES__H
