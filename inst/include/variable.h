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
        void init_proposal(const std::string &proposal, double scale);

        /// Enforce user choice for constant.
        void enforce_user_const() noexcept { is_user_const_ = true; }
        /// If applicable, initialize value from prior.
        /// Reset default scale for adaptive proposal rule.
        void restart(rnd_gen &prng);
        /// Draw the value from the distribution (used to generate observation
        /// nodes during simulation; unlike restart it ignores the constant flag).
        void draw_value(rnd_gen &prng);
        /// Optimize attribute (see Attribute::optimize) and
        /// distribution arguments (see Function::optimize).
        void optimize() override final;
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
        void recover() noexcept
        {
            value_ = old_;
            accept_rate_ -= adapt_smoothing_;
        }

        double get_log_prior() const;
        double get_prior_log_ratio() const;
        double get_proposal_log_ratio() const;
        double get_move_log_ratio() const
        {
            switch (ratio_mode_)
            {
            case Ratio::Compute:
                return get_prior_log_ratio() - get_proposal_log_ratio();
            default: // case Ratio::None:
                return 0;
            }
        }

        /// Reset the acceptance rate estimate.
        /// NOTE: initialised at the target to avoid first steps adaptation.
        void reset_accept_rate() noexcept { accept_rate_ = TARGET_ACCEPT_RATE; }
        /// Adapt the proposal rule from the acceptance estimate.
        void adapt_proposal(double ratio);
        /// Set EMA smoothing factor (0, 1). NOTE: shared across variables.
        static void set_adapt_smoothing(double a) noexcept { adapt_smoothing_ = a; }
        /// Get EMA smoothing factor.
        static double get_adapt_smoothing() noexcept { return adapt_smoothing_; }

#ifndef NDEBUG
        ////////////////////
        //////  SHOW  //////
        ////////////////////

        /// Get a description about the functional status.
        std::string str_summary() const override final;
#endif

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

    public:
        //////////////////////////
        //////  ADAPTATION  //////
        //////////////////////////

        /// Target Metropolis acceptance rate.
        /// Based on @cite roberts_weak_1997
        static constexpr double TARGET_ACCEPT_RATE = 0.44;

    private:
        /// Exponentially moving average (EMA) of the acceptance indicator.
        double accept_rate_ = TARGET_ACCEPT_RATE;
        /// EMA smoothing factor. NOTE: effective window is 1. / adapt_smoothing_.
        inline static double adapt_smoothing_ = 1.;

        ///////////////////////
        //////  HELPERS  //////
        ///////////////////////

        /// Assert that range is non-degenerate.
        void helper_assert_non_degenerate_range_();

        /// Initialize null proposal.
        void helper_init_null_proposal_();
        /// Initialize proposal of direct sample from prior/model.
        void helper_init_sample_proposal_();
        /// Initialize parametric Monte Carlo proposal.
        void helper_init_monte_carlo_proposal_(const std::string &proposal, double scale);
    };
}
