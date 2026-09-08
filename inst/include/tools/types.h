/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#include <map>

namespace housepi
{
    template <typename T>
    using dict = std::map<std::string, T, std::less<>>;
    using str_table = dict<std::vector<std::string>>;
    using dbl_table = dict<std::vector<double>>;

    /// Global individual index.
    using IndivIdx = std::uint32_t;
    /// Global household index.
    using HouseIdx = std::uint32_t;
    /// Relative member index within household.
    using MemberIdx = std::uint8_t;
    /// Individual status level.
    using StatusLvl = std::uint8_t;
    /// Storage slot a status level reads in a risk store.
    using StatusSlot = std::uint8_t;
    /// Offset of a status slot within one risk definition.
    using SlotOffset = std::uint16_t;
    /// Offset into a flat risk store.
    using RiskOffset = std::uint32_t;
    /// Fixed-size table addressed by a status level.
    template <typename T>
    using StatusTable = std::array<T, 1 + std::numeric_limits<StatusLvl>::max()>;
}
