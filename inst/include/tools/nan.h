/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>

namespace housepi::na
{
    constexpr double DOUBLE = std::numeric_limits<double>::quiet_NaN(); /// NA double.
    constexpr std::size_t INDEX = static_cast<std::size_t>(-1);         /// NA index.
    constexpr std::uint32_t UINDEX = static_cast<std::uint32_t>(-1);    /// NA narrow index.

    /// Whether double `x` is NA.
    inline bool is(double x) noexcept { return std::isnan(x); }
    /// Whether string `txt` is NA.
    bool is(std::string_view txt) noexcept;
    /// Whether index `i` is NA.
    constexpr bool is(std::size_t i) noexcept { return i == INDEX; }
}

namespace housepi::inf
{
    constexpr double DOUBLE = std::numeric_limits<double>::infinity(); /// Infinite double.
}
