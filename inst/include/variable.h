/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
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

        Variable();
        ~Variable();

        /// Copy constructor
        Variable(const Variable &other);
        /// Move constructor
        Variable(Variable &&);

        /// Copy assignment
        Variable &operator=(Variable other);

        void swap(Variable &other) noexcept;

        /// Sugar to call variable as a function if applicable, otherwise raise error.
        double operator()(double x) const;

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Restore empty. Used in multiple initializations.
        void clear() noexcept;

        /// Enforce user choice for constant.
        void enforce_user_const() noexcept { is_user_const_ = true; }

        /// Change value. Store old value.
        void set(double x)
        {
            old_ = value_;
            Attribute::set(x);
        }

        /// Recover old value.
        void recover() noexcept { value_ = old_; }

        /// Assign an attribute to function.
        void init_func_attr(std::string_view key, const Attribute &attr);
        /// Finalize reshaping of function.
        /// Drop function if it's constant when reshaping.
        /// NOTE: discretization of model has priority
        ///       e. g. user discr is FALSE, but discrete-defined model
        void init_func_reshaping(bool discr);

        /// Initialize pointer to distribution.
        void init_distr(std::string_view opt);
        /// Assign an attribute to distribution.
        void init_distr_attr(std::string_view key, const Attribute &attr);
        /// Finalize reshaping of distribution.
        /// NOTE: to be called necessarily.
        void init_distr_reshaping(bool discr);

        /// Restore user defined variable initialization.
        void user_reset(rnd_gen &prng);
        /// Optimizing pruning indirections and setting constants to numeric.
        void optimize() override final;
        /// Initialize pointer to proposer.
        /// NOTE: it must be assigned after model initialization.
        void init_proposal(const std::string &proposal, double proposal_scale);

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Wheter the variable is discrete.
        bool is_discrete() const noexcept override final;
        /// Wheter the variable has a prior or a model associated.
        bool has_model() const noexcept;
        /// Wheter the variable has an adaptable parametric proposal rule.
        bool is_adaptable() const noexcept;

        ///////////////////////////////////
        //////  METROPOLIS-HASTINGS  //////
        ///////////////////////////////////

        /// Set new proposal internally and store old value.
        void propose(rnd_gen &prng);
        bool is_valid_move() const;

        double get_distr_log_pdf() const;
        double get_distr_log_ratio() const;
        double get_propose_log_ratio() const;
        double get_move_log_ratio() const
        {
            switch (ratio_mode_)
            {
            case Ratio::Compute:
                return get_distr_log_ratio() - get_propose_log_ratio();
            default: // case Ratio::None:
                return 0;
            }
        }

        void update_accept_rate(bool is_accepted);
        void reset_accept_rate();
        void adapt_proposal(double ratio);

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
        std::unique_ptr<Distribution> pDistr_; /// prior\model distribution of the variable
        std::unique_ptr<Proposal> pProp_;      /// Monte Carlo proposer

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
        void helper_init_monte_carlo_proposal_(const std::string &proposal, double proposal_scale);
    };
}
