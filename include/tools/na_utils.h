/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 * 
 * Author: Francesco Camaglia
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 * See the COPYRIGHT file for details.
 * 
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the COPYING file for details.
 *
 * tools/na_utils.h
 *
 */

#ifndef HOUSEPI__TOOLS__NA_UTILS__H
#define HOUSEPI__TOOLS__NA_UTILS__H

#include <limits>
#include <cmath>
#include <string>

namespace hsp
{
    namespace na
    {
        /// NA value for double.
        constexpr double DOUBLE = std::numeric_limits<double>::quiet_NaN();
        /// NA value for string.
        constexpr const char *STRING = "NA";

        /// Internal check for double representing NA.
        inline bool is(double x) { return std::isnan(x); }

        /// Internal check for string representing NA.
        inline bool is(const std::string &txt)
        {
            return txt.empty() ||
                   txt == "<NA>" ||
                   txt == "NA" ||
                   txt == "NaN" ||
                   txt == "NAN" ||
                   txt == "nan";
        }

        inline bool is(const char *txt)
        {
            if (!txt)
                return true;
            return is(std::string(txt));
        }
    }
}

#endif // HOUSEPI__TOOLS__TYPES__H
