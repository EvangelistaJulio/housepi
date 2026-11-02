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

        Variable();           /// Constructor.
        ~Variable() override; /// Destructor.
        /// Copy constructor.
        /// * proposal -> not copied, the copy binds its own.
        /// @note Copy only before proposals are bound, as State allocation does.
        Variable(const Variable &other);
        Variable(Variable &&) noexcept;      /// Move constructor.
        Variable &operator=(Variable other); /// Copy assignment.

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
        /// * `Unconstrained` -> proposal log-ratio finite, infinite on a value rounded onto a bound;
        /// * otherwise -> always.
        bool is_valid_move() const
        {
            switch (move_mode_)
            {
            case Move::RandomWalk:
                return distr_->is_valid(value_);
            case Move::Unconstrained:
                return proposal_log_ratio() < inf::DOUBLE;
            default: // case Move::Forward, Move::Independence:
                return true;
            }
        }

        /// Restore previous value.
        void restore() noexcept { value_ = old_; }

        /// Whether the variable can join a block: it moves by a random walk in the unconstrained coordinate.
        bool can_join_block() const noexcept { return move_mode_ == Move::Unconstrained; }
        /// Jump by `dz` in the unconstrained coordinate, keeping the previous value.
        /// @return Proposal log-ratio of the move.
        /// @note Expects can_join_block().
        double block_jump(double dz)
        {
            RndWalkProposal &prop = unconstrained_();
            set(prop.jump(value_, dz));
            return prop.log_ratio(old_, value_);
        }
        /// Current value in the unconstrained coordinate.
        /// @note Expects can_join_block().
        double block_coord() const { return unconstrained_().to_coord(value_); }
        /// Default jump scale in the unconstrained coordinate.
        /// @note Expects can_join_block().
        double block_scale() const { return unconstrained_().default_scale(distr_.get()); }

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
            Forward,      /// Never moves, or draws from its own distribution: no ratio, always valid.
            Independence, /// Draws inside the domain: ratio computed, always valid.
            RandomWalk,   /// May move outside the domain: ratio computed, value checked.
            Unconstrained /// Moves in the unconstrained coordinate: ratio computed, proposal judges validity.
        } move_mode_ = Move::Forward;

        double old_ = na::DOUBLE;             /// Previous value.
        std::unique_ptr<Distribution> distr_; /// Distribution.
        std::unique_ptr<Proposal> prop_;      /// Proposal.

        //////  Helpers  //////

        /// Proposal as a random walk in the unconstrained coordinate.
        RndWalkProposal &unconstrained_() { return static_cast<RndWalkProposal &>(*prop_); }
        const RndWalkProposal &unconstrained_() const { return static_cast<const RndWalkProposal &>(*prop_); }

        /// Initialize null proposal.
        void helper_init_null_proposal_();
        /// Initialize independence proposal drawing directly from distribution.
        void helper_init_forward_proposal_();
        /// Initialize a random-walk proposal: rw, rw_log or rw_unif.
        void helper_init_random_walk_proposal_(const std::string &proposal, double scale);
        /// Initialize the default proposal, a random walk in the unconstrained coordinate.
        /// * discrete prior -> rw on its lattice.
        void helper_init_default_proposal_();
        /// Initialize an independence proposal drawing from a named distribution.
        void helper_init_independence_proposal_(const std::string &proposal, const dict<Attribute> &args);
    };
}
