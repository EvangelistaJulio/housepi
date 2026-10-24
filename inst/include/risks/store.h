/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "data/individual.h"
#include "data/status.h"
#include "function.h"

//////////////////////////////  risks::Store  //////////////////////////////

namespace housepi
{
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

namespace housepi::risks
{
    /// Risk attributes of every unique definition, in one flat store.
    ///
    /// * reading one crosses a single indirection;
    /// * each owner maps to the offset its block starts at;
    /// * owner: member for group risks, household for contact risks.
    struct Store
    {
    public:
        /// Whether the last build received no definition.
        bool is_empty() const noexcept { return flat_.empty(); }
        /// Finalize every attribute (see Attribute::finalize).
        void finalize();
        /// Households reading each parameter, by parameter key, in increasing order.
        const std::vector<std::vector<HouseIdx>> &houses_of_param() const noexcept { return houses_of_param_; }

    protected:
        std::vector<Attribute> flat_ = {};                        /// Risk attributes, one block per unique definition.
        std::vector<RiskOffset> base_ = {};                       /// Offset each owner's block starts at.
        std::vector<std::vector<HouseIdx>> houses_of_param_ = {}; /// Households reading each parameter.

        /// Record that household `H` reads parameter `p`.
        /// @note Households are visited in increasing order.
        void mark_reads_(size_t p, HouseIdx H)
        {
            if (p >= houses_of_param_.size())
                houses_of_param_.resize(p + 1);
            if (houses_of_param_[p].empty() || houses_of_param_[p].back() != H)
                houses_of_param_[p].push_back(H);
        }

        /// Drop the previous build and point every owner at the first block.
        void reset_store_(size_t n_owners);

        /// Attribute `offset` cells into the block that `owner` reads.
        const Attribute &cell_(size_t owner, size_t offset) const noexcept
        {
            return flat_[base_[owner] + offset];
        }
    };
}
