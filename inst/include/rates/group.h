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
#include "tools/keys.h"
#include "data/individual.h"
#include "rates/anchor.h"
#include "rates/store.h"

//////////////////////////////  rates::Group  //////////////////////////////

namespace housepi
{
    namespace data
    {
        struct Wide;
    }
    struct Parameter;
}

namespace housepi::rates
{
    /// Group rates of every member, keyed by the model-declared group-rate keys.
    ///
    /// * one block per unique definition, one slot per status a rule keys on;
    /// * a member drawing a new status reads the matching slot, no storage rewritten.
    struct Group : public KeyHolder, public Store
    {
    public:
        using KeyHolder::KeyHolder;

        //////  Setter  //////

        /// Build group rates from rule table `input`, one block per unique member definition.
        /// * empty `input` -> every declared rate reads the reference.
        /// @note Run it even on an empty `input`, so declared rates are initialized.
        /// @return Data columns the rules read, as equivalence-class covariates.
        std::vector<std::string> build(data::Wide &data_wide,
                                       const Parameter &param,
                                       const str_table &input);

        /// Declare the keys `vec` a rule may stratify by infected status.
        /// * a rate read on a member without infected status is excluded.
        void request_status_stratified(std::vector<size_t> vec) { status_stratified_ = std::move(vec); }

        //////  Getter  //////

        /// Keys a rule may stratify by infected status.
        const std::vector<size_t> &status_stratified() const noexcept { return status_stratified_; }

        /// Whether members `a` and `b` share a group-rate partition.
        /// * members reading different status slots read different rates.
        bool same_group(const Individual &a,
                        const Individual &b) const noexcept
        {
            return base_[Index::of(a)] == base_[Index::of(b)] &&
                   slot_base_[a.status()] == slot_base_[b.status()];
        }

        /// Group-rate attribute `g` of member `indiv`.
        /// * the status is folded into the address: a member drawing a new one reads the matching cell.
        const Attribute &at(size_t g, const Individual &indiv) const noexcept
        {
            return cell_(Index::of(indiv), slot_base_[indiv.status()] + g);
        }
        /// Group-rate value `g` of member `indiv`.
        double value(size_t g, const Individual &indiv) const
        {
            return at(g, indiv).value();
        }

        /// Rule groups of key `g`, their scale anchoring.
        const std::vector<ScaleAnchor> &scale_anchors(size_t g) const { return scale_anchors_[g]; }

        /// Number of rule groups of key `g` whose overall scale is free.
        size_t n_scale_free(size_t g, const Parameter &param) const
        {
            return count_scale_free(scale_anchors_[g], param);
        }

    private:
        std::vector<size_t> status_stratified_ = {}; /// Keys stratifiable by status.
        size_t n_slots_ = 1;                         /// Status slots stored per group rate.
        /// Offset each status level's slot starts at, pre-multiplied: a read costs one load and an add.
        /// * levels no rule distinguishes -> 0.
        StatusTable<SlotOffset> slot_base_ = {};
        std::vector<std::vector<ScaleAnchor>> scale_anchors_ = {}; /// Rule groups per group-rate key.

        /// Throw if a rule keys by status a rate the model does not allow.
        void assert_status_stratified_(const str_table &input,
                                       const std::vector<std::string> &g_input_vec) const;
        /// Give every status a rule keys on a slot of its own, leaving the rest on slot 0.
        /// @return Slot of each status level, as rule resolution reads it.
        StatusTable<StatusSlot> helper_init_status_slots_(const std::vector<std::string> &status_col);
        /// Register the rate names of `input`.
        /// @return Covariate columns a rule may key on.
        std::vector<std::string> helper_process_definitions_(const str_table &input);
        /// Append the block of one unique definition, a product of parameters per slot and key.
        /// @param cell_params Parameter keys per cell, slot-major.
        void helper_init_unique_rate_(const Parameter &param,
                                      const std::vector<std::vector<size_t>> &cell_params);
    };
}
