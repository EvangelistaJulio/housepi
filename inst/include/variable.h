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

#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Variable : class *** ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    class Distribution;
    class Proposal;

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
        void init_func_attr(std::string_view key, const Attribute &attr);
        /// Finalize reshaping of function.
        /// Drop function if it's constant when reshaping.
        /// NOTE: discretization of model has priority
        /// e. g. user discr is FALSE, but discrete-defined model
        void init_func_reshaping(bool discr);

        /// Initialize pointer to distribution.
        void init_prior(std::string_view opt);
        /// Assign an attribute to distribution.
        void init_prior_attr(std::string_view key, const Attribute &attr);
        /// Finalize reshaping of distribution.
        /// NOTE: to be called necessarily.
        void init_prior_reshaping(bool discr);
        /// Initialize pointer to proposal.
        /// NOTE: it must be assigned after model initialization.
        /// Random-walk kernels (rw/rw_log/rw_unif) read a "scale" argument; independence
        /// samplers (prior, or a distribution name) draw from an owned distribution.
        void init_proposal(const std::string &proposal, const dict<Attribute> &args);

        /// Enforce user choice for constant.
        void enforce_user_const() noexcept { is_user_const_ = true; }
        /// If applicable, initialize value from prior.
        /// Reset default scale for adaptive proposal rule.
        void restart(rnd_gen &prng);
        /// Draw the value from the distribution (used to generate observation
        /// nodes during simulation; unlike restart it ignores the constant flag).
        void draw_value(rnd_gen &prng);
        /// Finalize attribute (see Attribute::finalize) and
        /// distribution arguments (see Function::finalize).
        void finalize() override final;
        /// Restore empty. Used in multiple initializations.
        void clear() noexcept;

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Wheter the variable is discrete.
        bool is_discrete() const noexcept override final;
        /// Wheter the variable has a prior or a model associated.
        bool has_prior() const noexcept;
        /// Wheter the variable has an adaptive parametric proposal rule.
        bool is_adaptive() const noexcept;

        ///////////////////////////////////
        //////  METROPOLIS-HASTINGS  //////
        ///////////////////////////////////

        /// Set new proposal internally and store old value.
        /// NOTE: tentatively credits the move as accepted.
        void propose(rnd_gen &prng);
        /// Check whether last propose move stored a valid value.
        bool is_valid_move() const;
        /// Recover old value and roll back the tentative accept rate update.
        void recover() noexcept;

        /// Log-pdf of the inherent distribution field.
        /// It assumes current value field is valid (used to avoid double checks).
        double log_pdf() const;
        /// Guarded log-pdf guarded.
        /// Returns -Inf when the value is not valid (used when move not granted valid).
        double guarded_log_pdf() const;
        /// Log-ratio between new and old values.
        /// It assumes the pdf arguemnts (included support) did not change.
        double log_ratio() const;
        double proposal_log_ratio() const;
        /// Log-ratio of priors, minus log-ratio of proposals for last move.
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

        /// Reset the acceptance rate estimate.
        void reset_accept_rate() noexcept;
        /// Adapt the proposal rule from the acceptance estimate.
        void adapt_proposal(double ratio);
        static void set_adapt_smoothing(double s) noexcept;
        static double get_adapt_smoothing() noexcept;

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
        std::unique_ptr<Distribution> pDistr_; /// prior\model distribution
        std::unique_ptr<Proposal> pProp_;      /// Monte Carlo proposal

    private:
        ///////////////////////
        //////  HELPERS  //////
        ///////////////////////

        /// Assert that range is non-degenerate.
        void helper_assert_non_degenerate_range_();

        /// Initialize null proposal.
        void helper_init_null_proposal_();
        /// Initialize independence proposal drawing directly from the prior.
        void helper_init_prior_proposal_();
        /// Initialize a random-walk proposal (rw/rw_log/rw_unif kernel).
        void helper_init_random_walk_proposal_(const std::string &proposal, double scale);
        /// Initialize an independence proposal drawing from a named distribution.
        void helper_init_independence_proposal_(const std::string &proposal, const dict<Attribute> &args);
    };
}
