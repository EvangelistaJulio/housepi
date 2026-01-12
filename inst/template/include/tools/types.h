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

#include <stdexcept>
#include <string>
#include <vector>
#include <unordered_map>
#include <limits>

using int_vec = std::vector<int>;
using dbl_vec = std::vector<double>;
using str_vec = std::vector<std::string>;

/// Matrices
template <typename T>
using matrix = std::vector<std::vector<T>>;

/// Unordered dictionaries
template <typename T>
using u_dict = std::unordered_map<std::string, T>;

namespace hsp
{
    namespace na
    {
        constexpr double DOUBLE = std::numeric_limits<double>::quiet_NaN(); /// NA double value
        constexpr const char *STRING = "NA";                                /// NA string value

        /// Wheter double represents NA.
        bool is(double x);
        /// Wheter string represents NA.
        bool is(const std::string &txt);
        /// Wheter character represents NA.
        bool is(const char *txt);
    }
}

#endif // HOUSEPI__TOOLS__TYPES__H
