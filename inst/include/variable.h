/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "tools/dictionary.h"
#include "tools/nan.h"
#include "tools/random.h"
#include "distribution.h"
#include "function.h"
#include "proposal.h"

//////////////////////////////  Variable  //////////////////////////////

namespace housepi
{
    /// Attribute drawn from a distribution and moved by a proposal.
    class Variable final : public Attribute
    {
    public:
        //////  Constructor  //////

        Variable();
        ~Variable() override;
        /// Copy constructor.
        /// * proposal -> not copied, the copy binds its own.
        /// @note Copy only before proposals are bound, as State allocation does.
        Variable(const Variable &other);
        Variable(Variable &&) noexcept;
        Variable &operator=(Variable other);

        /// Swap contents with `other`.
        void swap(Variable &other) noexcept;

        //////  Setter  //////

        /// Set value to `x`, keeping the previous one.
        void set(double x)
        {
            old_ = value_;
            Attribute::set(x);
        }

        /// Reshape the function, see Function::reshape().
        /// @note The model's discretization takes priority, e.g. `discr` false on a discrete-defined model.
        void reshape_function(bool discr) { func_->reshape(discr); }

        /// Initialize distribution `opt`.
        void init_distr(std::string_view opt, bool discr);
        /// Initialize argument `key` of the distribution from `attr`.
        void init_distr_attr(std::string_view key, const Attribute &attr) { distr_->init_attr(key, attr); }
        /// Reshape the distribution.
        /// @note Always call it after init_distr().
        void reshape_distribution(bool discr);
        /// Initialize proposal `proposal` with arguments `args`.
        ///
        /// * random-walk kernels (rw, rw_log, rw_unif) -> read a `scale` argument;
        /// * independence samplers -> draw from an owned distribution.
        /// @note Call after the model is initialized.
        void init_proposal(const std::string &proposal, const dict<Attribute> &args);

        /// Enforce user choice for constant.
        void enforce_user_const() noexcept { is_user_const_ = true; }
        /// Drop user choice for constant. Used when readers bind its variable.
        void release_user_const() noexcept { is_user_const_ = false; }
        /// Pin within `range`: no density, no move, `range` reported as its values.
        void pin_to(Codomain range);
        /// Draw the initial value from the distribution, if any.
        void restart(RndGen &prng);
        /// Draw the value from the distribution, the constant flag ignored.
        /// * no distribution -> skipped.
        void draw_value(RndGen &prng)
        {
            if (distr_)
                set(distr_->rand(prng));
        }
        /// Finalize the attribute and the distribution arguments, see Attribute::finalize().
        void finalize() override final;
        /// Restore empty, between initializations.
        void clear() noexcept;

        //////  Flags  //////

        /// Values the variable takes, in the absolute frame.
        /// * constant -> its value, whatever the distribution;
        /// * otherwise -> the support of the distribution.
        Codomain codomain() const noexcept override final
        {
            return !distr_ || is_constant() ? Attribute::codomain() : distr_->support();
        }
        /// Whether the variable has a distribution.
        bool has_distr() const noexcept { return distr_ != nullptr; }
        /// Whether the variable has a definition: a value, a pointer, a function or a distribution.
        bool has_definition() const noexcept { return !is_empty() || distr_; }

        //////  Metropolis-Hastings  //////

        /// Propose a new value, keeping the previous one.
        void propose(RndGen &prng) { set(prop_->propose(prng, value_)); }
        /// Whether the last proposed value is valid.
        /// * `RandomWalk` -> inside the domain of the distribution;
        /// * `ProposalChecked` -> proposal log-ratio finite, infinite on values rounded onto a bound or off the domain;
        /// * otherwise -> always.
        bool is_valid_move() const
        {
            switch (move_mode_)
            {
            case Move::RandomWalk:
                return distr_->is_valid(value_);
            case Move::ProposalChecked:
                return proposal_log_ratio() < inf::DOUBLE;
            // case Move::Forward, Move::Independence:
            default:
                return true;
            }
        }

        /// Restore previous value.
        void restore() noexcept { value_ = old_; }
        /// Whether the last move proposed the value held, as a lattice jump rounded to 0.
        bool is_null_move() const noexcept { return value_ == old_; }

        /// Whether the variable can join a block: it moves by the default walk.
        bool can_join_block() const noexcept { return move_mode_ == Move::ProposalChecked; }
        /// Jump by `dz` in the coordinate of its walk, keeping the previous value.
        /// @return Proposal log-ratio of the move.
        /// @note Expects can_join_block().
        double block_jump(double dz)
        {
            RndWalkProposal &prop = walk_();
            set(prop.jump(value_, dz));
            return prop.log_ratio(old_, value_);
        }
        /// Current value in the coordinate of its walk.
        /// @note Expects can_join_block().
        double block_coord() const { return walk_().to_coord(value_); }
        /// Default jump scale in the coordinate of its walk.
        /// @note Expects can_join_block().
        double block_scale() const { return walk_().default_scale(distr_.get()); }
        /// Whether the jumps are rounded onto a unit lattice, see AdaptableProposal::is_lattice().
        /// @note Expects can_join_block().
        bool block_on_lattice() const noexcept { return walk_().is_lattice(); }

        /// Name of the proposal, see Proposal::name().
        std::string proposal_name() const { return prop_->name(); }
        /// Coordinate of the proposal jumps, see Proposal::coord().
        std::string proposal_coord() const { return prop_->coord(); }
        /// Current scale of the proposal, see Proposal::scale().
        double proposal_scale() const noexcept { return prop_->scale(); }

        /// Reset an adaptable proposal to its default scale.
        void init_proposal() { prop_->reset(distr_.get()); }
        /// Adapt the proposal to acceptance probability `alpha` of the last move, with gain `gain`.
        void adapt_proposal(double alpha, double gain) { prop_->adapt(alpha, gain); }

        /// Log-pdf of the current value.
        /// @note Assumes the value is valid: no check.
        double log_pdf() const { return distr_->log_pdf(value_); }
        /// Log-pdf of the current value, minus infinity when invalid.
        double guarded_log_pdf() const { return distr_->guarded_log_pdf(value_); }
        /// Log-ratio of densities from the previous to the current value.
        /// @note Assumes the density arguments, domain included, did not change.
        double log_ratio() const { return distr_->log_ratio(old_, value_); }
        /// Log-ratio of proposals from the previous to the current value.
        double proposal_log_ratio() const { return prop_->log_ratio(old_, value_); }
        /// Log-ratio of the last move: densities minus proposals.
        /// * `Forward` -> 0, the two cancel.
        double get_move_log_ratio() const
        {
            return move_mode_ == Move::Forward ? 0. : log_ratio() - proposal_log_ratio();
        }

    private:
        /// Proposal family, read by is_valid_move() and get_move_log_ratio().
        enum class Move : std::uint8_t
        {
            /// Never moves, or draws from its own distribution: no ratio, always valid.
            Forward,
            /// Draws inside the domain: ratio computed, always valid.
            Independence,
            /// May move outside the domain: ratio computed, value checked.
            RandomWalk,
            /// Default walk: ratio computed, proposal judges validity.
            ProposalChecked
        } move_mode_ = Move::Forward;

        /// Previous value.
        double old_ = na::DOUBLE;
        std::unique_ptr<Distribution> distr_;
        std::unique_ptr<Proposal> prop_;

        //////  Helpers  //////

        /// Proposal as the default random walk.
        RndWalkProposal &walk_() { return static_cast<RndWalkProposal &>(*prop_); }
        const RndWalkProposal &walk_() const { return static_cast<const RndWalkProposal &>(*prop_); }

        /// Initialize independence proposal drawing directly from distribution.
        void helper_init_forward_proposal_();
        /// Initialize a random-walk proposal: rw, rw_log or rw_unif.
        void helper_init_random_walk_proposal_(const std::string &proposal, double scale);
        /// Initialize the default proposal, a random walk in the coordinate mapping the domain onto the real line.
        /// * discrete prior -> walk on its lattice.
        void helper_init_default_proposal_();
        /// Initialize an independence proposal drawing from a named distribution.
        void helper_init_independence_proposal_(const std::string &proposal, const dict<Attribute> &args);
    };
}
