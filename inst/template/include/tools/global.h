/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__GLOBAL__H
#define HOUSEPI__TOOLS__GLOBAL__H

#include <string_view>

namespace hsp
{
    namespace discr
    {
        constexpr unsigned short NONE = 0;    /// None
        constexpr unsigned short LEFT = 1;    /// Left
        constexpr unsigned short CENTER = 2;  /// Center
        constexpr unsigned short DEFINED = 3; /// Defined
    }

    namespace infprof
    {
        /// Augmented instantaneous infectivity profile.
        constexpr const char *INSTANT = ".profile_inst";
        /// Augmented cumulative infectivity profile.
        constexpr const char *CUMUL = ".profile_cumul";
        /// Default size for lookup table for discrete case
        constexpr unsigned LOOKUP_SIZE = 30;
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

    namespace sources
    {
        /// Wheter definitions is 'augmented'.
        inline bool is_augmented(std::string_view s) { return s == "A"; }
        /// Wheter definitions is 'data'.
        inline bool is_data(std::string_view s) { return s == "D"; }
        /// Wheter definitions is 'group rate'.
        inline bool is_group_rate(std::string_view s) { return s == "G"; }
        /// Wheter definitions is 'numeric'.
        inline bool is_numeric(std::string_view s) { return s == "N"; }
        /// Wheter definitions is 'parameter'.
        inline bool is_parameter(std::string_view s) { return s == "P"; }
    }
}

#endif // HOUSEPI__TOOLS__GLOBAL__H
