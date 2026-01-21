/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

namespace hsp
{
    namespace discr
    {
        constexpr unsigned short NONE = 0;    /// None
        constexpr unsigned short LEFT = 1;    /// Left
        constexpr unsigned short CENTER = 2;  /// Center
        constexpr unsigned short DEFINED = 3; /// Defined
    }

    namespace params
    {
        constexpr const char *X = "x";     /// Evaluation parameter.
        constexpr const char *LOC = "loc"; /// Location parameter.
        constexpr const char *MIN = "min"; /// Minimum value.
        constexpr const char *MAX = "max"; /// Maximum value.
    }

    namespace rates
    {
        constexpr const char *REFERENCE = ".reference";  /// Rate reference choice.
        constexpr const char *MISSING = ".missing_rate"; /// Missing rate.
        constexpr const char *CONTACT = ".contact";      /// Contanct rate.
        constexpr double VALUE_REFERENCE = 1.0;          /// Alias for reference value.
    }
}
