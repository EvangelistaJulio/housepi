/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <limits>
#include <cstddef>

namespace housepi::na
{
    constexpr double DOUBLE = std::numeric_limits<double>::quiet_NaN(); /// NA double value
    constexpr std::size_t INDEX = static_cast<std::size_t>(-1);         /// NA index value

    /// Wheter double represents NA.
    bool is(double x) noexcept;
    /// Wheter string represents NA.
    bool is(std::string_view txt) noexcept;
    /// Wheter index represents NA.
    constexpr bool is(std::size_t i) noexcept { return i == INDEX; }
}
