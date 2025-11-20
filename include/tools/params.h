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
 * tools/params.h
 *
 */

#ifndef HOUSEPI__TOOLS__PARAMS__H
#define HOUSEPI__TOOLS__PARAMS__H

namespace hsp
{
    // Parameter names.
    namespace params
    {
        constexpr const char *X = "x";          /// Evaluation parameter.
        constexpr const char *LOCATION = "loc"; /// Location parameter.
        constexpr const char *MIN = "min";      /// Minimum value.
        constexpr const char *MAX = "max";      /// Maximum value.
    }
}

#endif // HOUSEPI__TOOLS__PARAMS__H
