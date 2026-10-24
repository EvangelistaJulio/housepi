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
#include "risks/anchor.h"
#include "risks/contact.h"
#include "risks/group.h"

//////////////////////////////  Risks  //////////////////////////////

namespace housepi
{
    namespace data
    {
        struct Wide;
    }
    struct Parameter;

    /// Group and contact risks of a run.
    ///
    /// * rules, resolved once per build into flat stores;
    /// * identifiability anchors, declared by the model.
    struct Risks
    {
        //////  Constructor  //////

        /// Register model-declared group-risk keys `keys`.
        explicit Risks(std::vector<std::string> keys)
            : group_(std::move(keys)) {}

        //////  Setter  //////

        /// Declare the group-risk keys `vec` a rule may stratify by infected status.
        void request_status_stratified(std::vector<size_t> vec) { group_.request_status_stratified(std::move(vec)); }
        /// Declare that baseline parameter `baseline` and group risks `keys` share one reference.
        /// @param contacts Whether contact risks scale the same term.
        void request_identifiability(size_t baseline, std::vector<size_t> keys, bool contacts)
        {
            anchors_.push_back({baseline, std::move(keys), contacts});
        }

        /// Build group risks from rule table `input` (see risks::Group::build).
        void build_group(data::Wide &data_wide, const Parameter &param, const str_table &input)
        {
            group_.build(data_wide, param, input);
        }
        /// Build contact risks from rule table `input` (see risks::Contact::build).
        void build_contact(data::Wide &data_wide, const Parameter &param, const str_table &input)
        {
            contact_.build(data_wide, param, input);
        }

        /// Finalize both stores (see Attribute::finalize).
        void finalize()
        {
            group_.finalize();
            contact_.finalize();
        }

        //////  Getter  //////

        /// Group risks.
        const risks::Group &group() const noexcept { return group_; }
        /// Contact risks.
        const risks::Contact &contact() const noexcept { return contact_; }
        /// Baselines the model anchors to risks.
        const std::vector<risks::Anchor> &anchors() const noexcept { return anchors_; }

        /// Whether members `i` and `j` of `house` are equivalent according to risks.
        /// * group risks match;
        /// * each contact risk of one matches its counterpart for the other, outgoing and incoming.
        bool same_partition(const Household &house, MemberIdx i, MemberIdx j) const;

        /// List the households reading each of `n_params` parameters through a risk, in increasing order.
        std::vector<std::vector<HouseIdx>> list_houses_of_param(size_t n_params) const;
        /// List the parameters moving together with each sampled baseline, the baseline first.
        /// * risk holding enough units against its reference -> joins the block of its baseline;
        /// * unknown units -> left out;
        /// * not sampled, or a declared proposal -> left out;
        /// * prior or bounds reading another member -> moves alone;
        /// * blocks sharing a member -> merged, a block of fewer than two members dropped.
        /// @note "Enough" is `risks::BLOCK_MIN_SHARE` of the units of the reference of its rule group.
        std::vector<std::vector<size_t>> list_param_blocks(const Parameter &param) const;
        /// Throw if more than one scale is free among a baseline and the risks scaling it.
        void assert_identifiability(const Parameter &param) const;

    private:
        risks::Group group_;                 /// Group risks, per member and status slot.
        risks::Contact contact_;             /// Contact risks, per ordered member pair.
        std::vector<risks::Anchor> anchors_; /// Baselines the model anchors to risks.
    };
}
