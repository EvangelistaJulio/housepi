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
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * tools/rates.h
 *
 */

#ifndef HOUSEPI__TOOLS__RATES__H
#define HOUSEPI__TOOLS__RATES__H

namespace hsp
{
    namespace rates
    {
        constexpr const char *REFERENCE = ".reference";  /// Rate reference choice.
        constexpr const char *MISSING = ".missing_rate"; /// Missing rate.
        constexpr const char *CONTACT = ".contact";      /// Contanct rate.
        constexpr double VALUE_REFERENCE = 1.0;          /// Alias for reference value.
    }
}

#endif // HOUSEPI__TOOLS__RATES__H
