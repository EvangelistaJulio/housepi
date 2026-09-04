/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/nan.h"
#include "tools/random.h"

#include "distribution.h"
#include "function.h"
#include "proposal.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Variable : class *** ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    class Variable final : public Attribute
    {
    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        Variable();                          /// constructor
        ~Variable();                         /// destructor
        Variable(const Variable &other);     /// copy
        Variable(Variable &&);               /// move
        Variable &operator=(Variable other); /// copy-assignment

        /// Swap method.
        void swap(Variable &other) noexcept;
        /// Sugar to call variable as a function if applicable, otherwise raise error.
        double operator()(double x) const;

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Change value. Store old value.
        void set(double x)
        {
            old_ = value_;
            Attribute::set(x);
        }

        /// Assign an attribute to function.
        void init_func_attr(std::string_view key, const Attribute &attr) { pFunc_->init_attr(key, attr); }
        /// Finalize reshaping of function.
        /// Drop function if it's constant when reshaping.
        /// NOTE: discretization of model has priority
        /// e. g. user discr is FALSE, but discrete-defined model
        void reshape_function(bool discr) { pFunc_->reshape(discr); }

        /// Initialize pointer to distribution.
        void init_distr(std::string_view opt, bool discr);
        /// Assign an attribute to distribution.
        void init_distr_attr(std::string_view key, const Attribute &attr) { pDistr_->init_attr(key, attr); }
        /// Finalize reshaping of distribution.
        /// NOTE: to be called necessarily.
        void reshape_distribution(bool discr);
        /// Initialize pointer to proposal.
        /// NOTE: it must be assigned after model initialization.
        /// Random-walk kernels (rw/rw_log/rw_unif) read a "scale" argument;
        /// independence samplers draw from an owned distribution.
        void init_proposal(const std::string &proposal, const dict<Attribute> &args);

        /// Enforce user choice for constant.
        void enforce_user_const() noexcept { is_user_const_ = true; }
        /// If applicable, initialize value from distribution.
        void restart(rnd_gen &prng);
        /// Draw the value from the distribution. Ignores the constant flag.
        /// A variable pinned by the simulation-only override is skipped.
        void draw_value(rnd_gen &prng)
        {
            // NOTE: distribution is dropped if pinned by the simulation
            if (has_distr())
                set(pDistr_->rand(prng));
        }
        /// Finalize attribute (see Attribute::finalize) and
        /// distribution arguments (see Function::finalize).
        void finalize() override final;
        /// Restore empty. Used in multiple initializations.
        void clear() noexcept;

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Values the variable takes, in the absolute frame.
        Codomain codomain() const noexcept override final
        {
            return has_distr() ? pDistr_->support() : Attribute::codomain();
        }
        /// Wheter the variable has a distribution associated.
        bool has_distr() const noexcept { return pDistr_ != nullptr; }
        /// Wheter the variable has an adaptive parametric proposal rule.
        bool is_adaptive() const noexcept { return pProp_ && pProp_->is_adaptive(); }

        ///////////////////////////////////
        //////  METROPOLIS-HASTINGS  //////
        ///////////////////////////////////

        /// Set new proposal internally and store old value.
        void propose(rnd_gen &prng) { set(pProp_->propose(prng, value_)); }
        /// Check whether last propose move stored a valid value.
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

        /// Recover old value.
        void restore() noexcept { value_ = old_; }

        /// Eventually reset the adaptable proposal to its default scale.
        void init_proposal() { pProp_->reset(pDistr_.get()); }
        /// Eventually adapt the proposal rule.
        void adapt_proposal(double ratio) { pProp_->adapt(ratio); }

        /// Setter/getter of shared smoothing factor for moving average in adaptable proposals.
        static void set_adapt_smoothing(double s) noexcept { AdaptableProposal::set_smoothing(s); }
        static double get_adapt_smoothing() noexcept { return AdaptableProposal::get_smoothing(); }

        /// Under adaptation, inform the proposal about refusal.
        void refuse_proposal_move() noexcept { pProp_->refuse_move(); }
        /// Under adaptation, inform the proposal about acceptance.
        void commit_proposal_move() noexcept { pProp_->commit_move(old_, value_); }

        /// Log-pdf of the inherent distribution field.
        /// It assumes current value field is valid (used to avoid double checks).
        double log_pdf() const { return pDistr_->log_pdf(value_); }
        /// Guarded log-pdf guarded.
        /// Returns -Inf when the value is not valid (used when move not granted valid).
        double guarded_log_pdf() const { return pDistr_->guarded_log_pdf(value_); }
        /// Log-ratio between new and old values.
        /// It assumes the pdf arguments (domain included) did not change.
        double log_ratio() const { return pDistr_->log_ratio(old_, value_); }
        double proposal_log_ratio() const { return pProp_->log_ratio(old_, value_); }
        /// Log-ratio of distributions, minus log-ratio of proposals for last move.
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
        enum class Validity : std::uint8_t
        {
            Always,
            Check
        } valid_mode_ = Validity::Always;
        enum class Ratio : std::uint8_t
        {
            None,
            Compute
        } ratio_mode_ = Ratio::None;

        double old_ = na::DOUBLE;              /// old value
        std::unique_ptr<Distribution> pDistr_; /// defining distribution
        std::unique_ptr<Proposal> pProp_;      /// Monte Carlo proposal

    private:
        ///////////////////////
        //////  HELPERS  //////
        ///////////////////////

        /// Initialize null proposal.
        void helper_init_null_proposal_();
        /// Initialize independence proposal drawing directly from distribution.
        void helper_init_distr_proposal_();
        /// Initialize a random-walk proposal (rw/rw_log/rw_unif kernel).
        void helper_init_random_walk_proposal_(const std::string &proposal, double scale);
        /// Initialize an independence proposal drawing from a named distribution.
        void helper_init_independence_proposal_(const std::string &proposal, const dict<Attribute> &args);
    };
}
