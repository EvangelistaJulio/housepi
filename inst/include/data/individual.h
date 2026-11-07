/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>

#include "data/status.h"

//////////////////////////////  Individual  //////////////////////////////

namespace housepi
{
    /// Global individual index.
    using IndivIdx = std::uint32_t;
    /// Global household index.
    using HouseIdx = std::uint32_t;
    /// Household equivalence class index.
    using EqClassIdx = std::uint32_t;
    /// Member index within its household.
    using MemberIdx = std::uint8_t;

    /// Member of one household.
    ///
    /// The user id lives in data::Dataset, keyed by individual index, see Index::of().
    class Individual final
    {
        /// Owns and mutates its members.
        friend class Household;
        /// Reads the opaque storage indexes.
        friend struct Index;

    public:
        //////  Operator  //////

        /// Whether they are different individuals.
        bool operator!=(const Individual &other) const noexcept { return indiv_index_ != other.indiv_index_; }

        //////  Getter  //////

        /// Status.
        StatusLvl status() const noexcept { return status_; }

        /// Whether status is infected.
        bool is_infected_status() const noexcept { return status::is_infected(status_); }
        /// Whether status is missing.
        bool is_missing_status() const noexcept { return status::is_missing(status_); }

    private:
        /// Global individual index.
        IndivIdx indiv_index_ = 0;
        /// Index within its household.
        MemberIdx member_index_ = 0;
        StatusLvl status_ = status::MISSING;
    };

    //////////////////////////////  Index  //////////////////////////////

    class Household;

    /// Storage indexes of Individual and Household, a client-attorney for lightweight viewers.
    struct Index
    {
        /// Index of `indiv` within its household.
        static MemberIdx member(const Individual &indiv) noexcept { return indiv.member_index_; }
        /// Global index of `indiv`.
        static IndivIdx of(const Individual &indiv) noexcept { return indiv.indiv_index_; }
        /// Global index of `house`, defined in data/household.h.
        static HouseIdx of(const Household &house) noexcept;
    };
}
