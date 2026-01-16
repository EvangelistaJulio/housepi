/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__NAN__H
#define HOUSEPI__TOOLS__NAN__H

/// DEVELOPER: flag used to develop Rcpp chunks
// #define _HOUSEPI_RCPP

#include <string>
#include <limits>

namespace hsp::na
{
    constexpr double DOUBLE = std::numeric_limits<double>::quiet_NaN(); /// NA double value

    /// Wheter double represents NA.
    bool is(double x);
    /// Wheter string represents NA.
    bool is(std::string_view txt);
}

#endif // HOUSEPI__TOOLS__NAN__H
