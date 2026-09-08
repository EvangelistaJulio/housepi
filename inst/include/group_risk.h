/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "tools/keys.h"
#include "tools/types.h"

#include "function.h"
#include "individual.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** GroupRisk : struct ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct DataWide;
    struct Parameter;

    /// Scale anchoring of a single rule group (see Parameter::assert_risk_identifiability).
    /// Resolving the group against a member, or an ordered pair for contact risks, either
    /// matches a rule, and reads the risk it names, or matches none, and reads the reference.
    /// The risks read are kept rather than counted: whether they can move is answered
    /// against the pinning of the run being asserted.
    struct ScaleAnchor
    {
        bool anchored = false;                   /// some resolution reads the reference
        std::vector<size_t> scaling_params = {}; /// distinct risks the group applies

        /// Register the risk read by one resolution. An empty name means no rule matched.
        void update(std::string_view pstr, const Parameter &param);
        /// Whether the overall scale of the rule group is left free:
        /// no resolution reads the reference and every risk applied can move.
        /// A single fixed risk pins the level of the whole group.
        bool is_free(const Parameter &param) const;
    };

    /// Risk attributes of every unique definition, in one flat store: reading one
    /// crosses a single indirection. An owner - a member for group risks, a household
    /// for contact risks - maps to the offset the block it reads starts at.
    struct RiskStorage
    {
    public:
        /// Whether no definition was provided in the last build.
        bool is_empty() const noexcept { return flat_.empty(); }
        /// Finalize all attributes (see Attribute::finalize).
        void finalize();

    protected:
        /// risk attributes, one block per unique definition
        std::vector<Attribute> flat_ = {};
        /// offset the block of each owner starts at
        std::vector<RiskOffset> base_ = {};

        /// Drop the previous build and address every owner to the first block.
        void reset_store_(size_t n_owners);

        /// Attribute 'offset' cells into the block owner 'owner' reads.
        const Attribute &cell_(size_t owner, size_t offset) const noexcept
        {
            return flat_[base_[owner] + offset];
        }

        /// Number of rule groups whose overall scale is left free.
        static size_t count_scale_free_(const std::vector<ScaleAnchor> &anchors,
                                        const Parameter &param);
    };

    struct GroupRisk : public KeyModule, public RiskStorage
    {
    public:
        /// Update model report and set group risks associate to new individual.
        /// NOTE: this has to run even if input is empty to initialize existing risks to reference.
        void build(DataWide &data_wide,
                   const Parameter &param,
                   const str_table &deftable);

        /// Declare the keys a rule may stratify by the infected status.
        /// A risk read on a member that carries no infected status is excluded.
        void request_status_stratified(std::vector<size_t> vec) { status_stratified_ = std::move(vec); }

        /// Whether two members belong to the same group-risk partition.
        /// Members reading different status slots read different risks.
        bool same_group(const Individual &a, const Individual &b) const noexcept
        {
            return base_[Index::of(a)] == base_[Index::of(b)] &&
                   slot_base_[a.status()] == slot_base_[b.status()];
        }

        /// Access the group risk attribute of given key for an individual.
        /// The status is folded into the address, so a member drawing a new one
        /// reads the matching cell with no storage rewritten.
        const Attribute &at(size_t g, const Individual &indiv) const noexcept
        {
            return cell_(Index::of(indiv), slot_base_[indiv.status()] + g);
        }
        /// Access the group risk value of given key for an individual.
        double value(size_t g, const Individual &indiv) const
        {
            return at(g, indiv).value();
        }

        /// Rule groups of given key whose overall scale is free (see Parameter::assert_risk_identifiability).
        size_t n_scale_free(size_t g, const Parameter &param) const
        {
            return count_scale_free_(scale_anchors_[g], param);
        }

    private:
        /// keys stratifiable by status
        std::vector<size_t> status_stratified_ = {};
        /// status slots stored per group risk
        size_t n_slots_ = 1;
        /// Status level to the offset its slot starts at, pre-multiplied so the
        /// read costs one load and an add. Undistinguished levels sit at 0.
        StatusTable<SlotOffset> slot_base_ = {};
        /// rule groups per group risk
        std::vector<std::vector<ScaleAnchor>> scale_anchors_ = {};

        /// Verify that only model allowed rules declare a key by status.
        void assert_status_stratified_(const str_table &input, const std::vector<std::string> &g_input_vec) const;
        /// Give every status a rule keys on its own slot, leaving the rest on slot 0.
        /// Returns the slot of each status level, as the rule resolution reads it.
        StatusTable<StatusSlot> helper_init_status_slots_(const std::vector<std::string> &status_col);
        /// Verify that group risks definitions are correct and pre-process theme.
        std::vector<std::string> helper_process_definitions_(DataWide &data_wide, const str_table &input);
        /// Convert to a dictionary of parameter pointers and pass the pointers to the individual.
        void helper_init_unique_risk_(const Parameter &param, const std::vector<str_table> &slot_risks);
    };
}
