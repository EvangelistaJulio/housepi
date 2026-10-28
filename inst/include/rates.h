/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "tools/dictionary.h"
#include "data/household.h"
#include "data/individual.h"
#include "rates/anchor.h"
#include "rates/contact.h"
#include "rates/group.h"

//////////////////////////////  Rates  //////////////////////////////

namespace housepi
{
    namespace data
    {
        struct Wide;
    }
    struct Parameter;

    /// Group and contact rates of a run.
    ///
    /// * rules, resolved once per build into flat stores;
    /// * identifiability anchors, declared by the model.
    struct Rates
    {
        //////  Constructor  //////

        /// Register model-declared group-rate keys `keys`.
        explicit Rates(std::vector<std::string> keys)
            : group_(std::move(keys)) {}

        //////  Setter  //////

        /// Declare the group-rate keys `vec` a rule may stratify by infected status.
        void request_status_stratified(std::vector<size_t> vec) { group_.request_status_stratified(std::move(vec)); }
        /// Declare that baseline parameter `baseline` and group rates `keys` share one reference.
        /// @param contacts Whether contact rates scale the same term.
        void request_identifiability(size_t baseline, std::vector<size_t> keys, bool contacts)
        {
            anchors_.push_back({baseline, std::move(keys), contacts});
        }

        /// Build group rates from rule table `input` (see rates::Group::build).
        std::vector<std::string> build_group(data::Wide &data_wide, const Parameter &param, const str_table &input)
        {
            return group_.build(data_wide, param, input);
        }
        /// Build contact rates from rule table `input` (see rates::Contact::build).
        std::vector<std::string> build_contact(data::Wide &data_wide, const Parameter &param, const str_table &input)
        {
            return contact_.build(data_wide, param, input);
        }

        /// Finalize both stores (see Attribute::finalize).
        void finalize()
        {
            group_.finalize();
            contact_.finalize();
        }

        //////  Getter  //////

        /// Group rates.
        const rates::Group &group() const noexcept { return group_; }
        /// Contact rates.
        const rates::Contact &contact() const noexcept { return contact_; }
        /// Baselines the model anchors to rates.
        const std::vector<rates::Anchor> &anchors() const noexcept { return anchors_; }

        /// Whether members `i` and `j` of `house` are equivalent according to rates.
        /// * group rates match;
        /// * each contact rate of one matches its counterpart for the other, outgoing and incoming.
        bool same_partition(const Household &house, MemberIdx i, MemberIdx j) const;

        /// List the households reading each of `n_params` parameters through a rate, in increasing order.
        std::vector<std::vector<HouseIdx>> list_houses_of_param(size_t n_params) const;
        /// List the parameters moving together with each sampled baseline, the baseline first.
        /// * rate holding enough units against its reference -> joins the block of its baseline;
        /// * unknown units -> left out;
        /// * not sampled, or a declared proposal -> left out;
        /// * prior or bounds reading another member -> moves alone;
        /// * blocks sharing a member -> merged, a block of fewer than two members dropped.
        /// @note "Enough" is `rates::BLOCK_MIN_SHARE` of the units of the reference of its rule group.
        std::vector<std::vector<size_t>> list_param_blocks(const Parameter &param) const;
        /// Throw if more than one scale is free among a baseline and the rates scaling it.
        void assert_identifiability(const Parameter &param) const;

    private:
        rates::Group group_;                 /// Group rates, per member and status slot.
        rates::Contact contact_;             /// Contact rates, per ordered member pair.
        std::vector<rates::Anchor> anchors_; /// Baselines the model anchors to rates.
    };
}
