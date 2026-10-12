/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
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

        Variable();                          /// Constructor.
        ~Variable() override;                /// Destructor.
        /// Copy constructor.
        /// - proposal: not copied, the copy binds its own.
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
        void reshape_function(bool discr) { pFunc_->reshape(discr); }

        /// Initialize distribution `opt`.
        void init_distr(std::string_view opt, bool discr);
        /// Initialize argument `key` of the distribution from `attr`.
        void init_distr_attr(std::string_view key, const Attribute &attr) { pDistr_->init_attr(key, attr); }
        /// Reshape the distribution.
        /// @note Always call it after init_distr().
        void reshape_distribution(bool discr);
        /// Initialize proposal `proposal` with arguments `args`.
        ///
        /// - random-walk kernels (rw, rw_log, rw_unif): read a `scale` argument;
        /// - independence samplers: draw from an owned distribution.
        /// @note Call after the model is initialized.
        void init_proposal(const std::string &proposal, const dict<Attribute> &args);

        /// Enforce user choice for constant.
        void enforce_user_const() noexcept { is_user_const_ = true; }
        /// Drop user choice for constant. Used when readers bind its variable.
        void release_user_const() noexcept { is_user_const_ = false; }
        /// Draw the initial value from the distribution, if any.
        void restart(rnd_gen &prng);
        /// Draw the value from the distribution, the constant flag ignored.
        /// - pinned by the simulation-only override: skipped.
        void draw_value(rnd_gen &prng)
        {
            // a variable pinned by the simulation holds no distribution
            if (pDistr_)
                set(pDistr_->rand(prng));
        }
        /// Finalize the attribute and the distribution arguments, see Attribute::finalize().
        void finalize() override final;
        /// Restore empty, between initializations.
        void clear() noexcept;

        //////  Flags  //////

        /// Values the variable takes, in the absolute frame.
        /// - constant: its value, whatever the distribution;
        /// - otherwise: the support of the distribution.
        Codomain codomain() const noexcept override final
        {
            return !pDistr_ || is_constant() ? Attribute::codomain() : pDistr_->support();
        }
        /// Whether the variable has a distribution.
        bool has_distr() const noexcept { return pDistr_ != nullptr; }
        /// Whether the variable has a definition: a value, a pointer, a function or a distribution.
        bool has_definition() const noexcept { return !is_empty() || pDistr_; }

        //////  Metropolis-Hastings  //////

        /// Propose a new value, keeping the previous one.
        void propose(rnd_gen &prng) { set(pProp_->propose(prng, value_)); }
        /// Whether the last proposed value is valid.
        bool is_valid_move() const
        {
            switch (valid_mode_)
            {
            case Validity::Check:
                return pDistr_->is_valid(value_);
            default: // case Validity::Always:
                return true;
            }
        }

        /// Restore previous value.
        void restore() noexcept { value_ = old_; }

        /// Reset an adaptable proposal to its default scale.
        void init_proposal() { pProp_->reset(pDistr_.get()); }
        /// Adapt the proposal to acceptance ratio `ratio`, if adaptable.
        void adapt_proposal(double ratio) { pProp_->adapt(ratio); }

        /// Under adaptation, inform the proposal about refusal, with accept-rate smoothing `s`.
        void refuse_proposal_move(double s) noexcept { pProp_->refuse_move(s); }
        /// Under adaptation, inform the proposal about acceptance, with accept-rate smoothing `s`.
        void commit_proposal_move(double s) noexcept { pProp_->commit_move(old_, value_, s); }

        /// Log-pdf of the current value.
        /// @note Assumes the value is valid: no check.
        double log_pdf() const { return pDistr_->log_pdf(value_); }
        /// Log-pdf of the current value, minus infinity when invalid.
        double guarded_log_pdf() const { return pDistr_->guarded_log_pdf(value_); }
        /// Log-ratio of densities from the previous to the current value.
        /// @note Assumes the density arguments, domain included, did not change.
        double log_ratio() const { return pDistr_->log_ratio(old_, value_); }
        /// Log-ratio of proposals from the previous to the current value.
        double proposal_log_ratio() const { return pProp_->log_ratio(old_, value_); }
        /// Log-ratio of the last move: densities minus proposals.
        double get_move_log_ratio() const
        {
            switch (ratio_mode_)
            {
            case Ratio::Compute:
                return log_ratio() - proposal_log_ratio();
            default: // case Ratio::None:
                return 0;
            }
        }

    private:
        /// Whether is_valid_move() checks the value.
        enum class Validity : std::uint8_t
        {
            Always,
            Check
        } valid_mode_ = Validity::Always;
        /// Whether get_move_log_ratio() computes the ratio.
        enum class Ratio : std::uint8_t
        {
            None,
            Compute
        } ratio_mode_ = Ratio::None;

        double old_ = na::DOUBLE;              /// Previous value.
        std::unique_ptr<Distribution> pDistr_; /// Distribution.
        std::unique_ptr<Proposal> pProp_;      /// Proposal.

    private:
        //////  Helpers  //////

        /// Initialize null proposal.
        void helper_init_null_proposal_();
        /// Initialize independence proposal drawing directly from distribution.
        void helper_init_distr_proposal_();
        /// Initialize a random-walk proposal: rw, rw_log or rw_unif.
        void helper_init_random_walk_proposal_(const std::string &proposal, double scale);
        /// Initialize an independence proposal drawing from a named distribution.
        void helper_init_independence_proposal_(const std::string &proposal, const dict<Attribute> &args);
    };
}
