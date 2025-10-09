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

        /// Internal check for std::string representing NA.
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
