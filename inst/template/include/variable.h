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
public:
    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    Variable(const std::string &err_type = "variable");
    ~Variable();

    /////////////////////////
    //////  OPERATORS  //////
    /////////////////////////

    /// Sugar to call variable as a function if applicable, otherwise raise error.
    double operator()(double x) const { return (this->*evaluate_)(x); }

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// Helper to initialize variable.
    void init(double value, bool is_const, std::string_view expr, bool is_prim = false);

    /// Recover buffered value.
    void recover_value() { value_ = buffer_; }

    /// Restore user defined variable initialization.
    void user_reset(hsp::rnd_gen &prng);

    /// Assign an attribute to function.
    void init_func_attr(std::string_view key, const Attribute &attr);

    /// Initialize function reshaping.
    /// NOTE: to be called necessarily.
    void init_func_reshaping(bool discr);

    /// Initialize pointer to distribution.
    void init_distr(std::string_view opt);

    /// Assign an attribute to distribution.
    void init_distr_attr(std::string_view key, const Attribute &attr);

    /// Finalize reshaping of distribution.
    /// NOTE: to be called necessarily.
    void init_distr_reshaping(bool discr);

    /// Random initialization of variable.
    void init_rand(hsp::rnd_gen &prng);

    /// Initialize pointer to proposer.
    void init_proposer(const std::string &proposal, double proposal_scale);

    //////////////////////
    //////  GETTER  //////
    //////////////////////

    /// Get pointer to distribution.
    const std::shared_ptr<Distribution> get_distribution() const { return pDistr_; }

    bool is_constant() const noexcept override final;
    bool is_discrete() const noexcept override final;
    bool is_model() const noexcept;
    bool is_adaptable() const noexcept;

    ///////////////////////////////////
    //////  METROPOLIS-HASTINGS  //////
    ///////////////////////////////////

    /// Set new proposal internally and store a buffer of old value.
    void propose(hsp::rnd_gen &prng);
    bool is_valid_move() const { return (this->*is_valid_move_)(); }

    double get_distr_logpdf() const;
    double get_distr_logratio() const;
    double get_propose_logratio() const;
    double get_move_logratio() const { return (this->*get_move_logratio_)(); }

    void update_accept_rate(bool is_accepted);
    void reset_accept_rate();
    void adapt_proposal(double ratio);

    ////////////////////
    //////  SHOW  //////
    ////////////////////

    /// Get a description about the functional status.
    std::string str_summary() const override final;

protected:
    double buffer_;           /// buffered value
    double user_value_;       /// initial user defined value
    bool is_constant_ = true; /// user flag for wheter constant or not

    // NOTE: Household must be copyable for simulations -> Variable pointers must be shared
    std::shared_ptr<Distribution> pDistr_; /// prior\model distribution of the variable
    std::shared_ptr<Proposal> pProp_;      /// Monte Carlo proposer

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Variable : getter [private]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double (Variable::*evaluate_)(double) const = &Variable::do_eval_error_;
    double do_eval_error_(double) const { throw std::runtime_error(
        "Called undefined evaluate operator for non-function " + err_name() + "."); }
    double do_eval_func_(double x) const;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Variable : setter [private]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Helper with proposer initialization of direct sample from prior/model.
    void helper_init_sample_proposer_();
    /// Helper with Monte Carlo proposer initialization.
    void helper_init_proposer_monte_carlo_(const std::string &proposal, double proposal_scale);

    /// Helper with degenerate range.
    void helper_degenerate_range_();

    /// Store the buffer and change value.
    void do_set_default_(double x) override final
    {
        buffer_ = value_;
        value_ = x;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  Variable : Metropolis-Hastings [private]  ///////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    //----- Valid move : functional methods -----//
    bool (Variable::*is_valid_move_)() const = &Variable::do_is_valid_move_default_;
    bool do_is_valid_move_default_() const;
    bool do_is_valid_move_always_() const { return true; }

    //----- Move logratio : functional methods -----//
    double (Variable::*get_move_logratio_)() const = &Variable::do_get_move_logratio_default_;
    double do_get_move_logratio_default_() const { return get_distr_logratio() - get_propose_logratio(); }
    double do_move_logratio_sample_() const { return 0.; }
};
