/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <map>

namespace housepi
{
    template <typename T>
    using dict = std::map<std::string, T, std::less<>>;
    using str_table = dict<std::vector<std::string>>;
    using dbl_table = dict<std::vector<double>>;

    using IndivIndex = std::uint32_t;  /// global individual index
    using HouseIndex = std::uint32_t;  /// global household index
    using MemberIndex = std::uint8_t;  /// index of a member within its household
}
