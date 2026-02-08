/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "tools/nan.h"
#include "tools/random.h"

#include "attribute.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Variable : class *** ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Distribution;
class Proposal;

class Variable : public Attribute
{
    enum class Validity : uint8_t
    {
        Always,
        Check
    } valid_mode_ = Validity::Check;
    enum class Ratio : uint8_t
    {
        None,
        Compute
    } ratio_mode_ = Ratio::Compute;

public:
    ~Variable();

    /// Sugar to call variable as a function if applicable, otherwise raise error.
    double operator()(double x) const;

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// Restore empty. Used in multiple initializations.
    void clear() noexcept
    {
        mode_ = Mode::Empty;
        value_ = hsp::na::DOUBLE;
        pAttr_ = nullptr;
        pFunc_ = nullptr;
        pDistr_ = nullptr;
        pProp_ = nullptr;
    }

    /// Enforce user choice for constant.
    void enforce_user_const() noexcept { is_user_const_ = true; }

    /// Change value. Buffer old value.
    void set(double x)
    {
        buffer_ = value_;
        Attribute::set(x);
    }

    /// Recover buffered value.
    void recover() noexcept { value_ = buffer_; }

    /// Initialize pointer to distribution.
    void init_distr(std::string_view opt);
    /// Assign an attribute to distribution.
    void init_distr_attr(std::string_view key, const Attribute &attr);
    /// Finalize reshaping of distribution.
    /// NOTE: to be called necessarily.
    void init_distr_reshaping(bool discr);

    /// Restore user defined variable initialization.
    void user_reset(hsp::rnd_gen &prng);
    /// Eventually prune indirections if to constants.
    void optimize();
    /// Initialize pointer to proposer.
    void init_proposal(const std::string &proposal, double proposal_scale);

    /// Assign an attribute to function.
    void init_func_attr(std::string_view key, const Attribute &attr);
    /// Initialize function reshaping.
    /// NOTE: to be called necessarily.
    void init_func_reshaping(bool discr);

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    bool is_discrete() const noexcept override final;
    bool is_model() const noexcept;
    bool is_adaptable() const noexcept;

    ///////////////////////////////////
    //////  METROPOLIS-HASTINGS  //////
    ///////////////////////////////////

    /// Set new proposal internally and store a buffer of old value.
    void propose(hsp::rnd_gen &prng);
    bool is_valid_move() const;

    double get_distr_logpdf() const;
    double get_distr_logratio() const;
    double get_propose_logratio() const;
    double get_move_logratio() const
    {
        switch (ratio_mode_)
        {
        case Ratio::None:
            return 0;
        case Ratio::Compute:
            return get_distr_logratio() - get_propose_logratio();
        }
        __builtin_unreachable();
    }

    void update_accept_rate(bool is_accepted);
    void reset_accept_rate();
    void adapt_proposal(double ratio);

    ////////////////////
    //////  SHOW  //////
    ////////////////////

    /// Get a description about the functional status.
    std::string str_summary() const override final;

protected:
    double buffer_ = hsp::na::DOUBLE;     /// buffered value

    /// NOTE: Household must be copyable for simulations 
    /// UPDATE: pointers should become unique -> Function clonable.
    std::shared_ptr<Distribution> pDistr_; /// prior\model distribution of the variable
    std::shared_ptr<Proposal> pProp_;      /// Monte Carlo proposer

    ///////////////////////
    //////  HELPERS  //////
    ///////////////////////

    /// Assert that range is non-degenerate.
    void helper_assert_non_degenerate_range_();
    /// Initialize proposal of direct sample from prior/model.
    void helper_init_sample_proposal_();
    /// Initialize parametric Monte Carlo proposal.
    void helper_init_monte_carlo_proposal_(const std::string &proposal, double proposal_scale);
};
