/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__VARIABLE__H
#define HOUSEPI__VARIABLE__H

#include <memory>

#include "tools/types.h"
#include "tools/na_utils.h"
#include "tools/random.h"
#include "tools/parser.h"
#include "function.h"
#include "distribution.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Proposal : class *** ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Base class for Metropolis proposers.
class Proposal
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Proposal : constructor  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Empty proposer.
    Proposal() {};

    /// Proposal constructor from Distribution pointer.
    Proposal(std::shared_ptr<Distribution>,
             double proposal_scale = hsp::na::DOUBLE,
             const std::vector<unsigned short> &args_to_adapt = {},
             bool is_discrete = false);

    /// Proposal constructor from Distribution name.
    Proposal(const std::string &distr_name,
             double proposal_scale = hsp::na::DOUBLE,
             const std::vector<unsigned short> &args_to_adapt = {},
             bool is_discrete = false);

    virtual ~Proposal() {};

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Proposal : getter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string name() const { return this->name_; }
    double get_arg(unsigned short argID) const { return pDistr_->arg(argID); }
    virtual double propose(rnd_gen &, double) { throw std::invalid_argument(
        "Called propose method of undefined Monte Carlo proposer."); }
    virtual double logratio(double, double) const { return 0.; }

    /// Wheter proposer can be adapted.
    bool is_adaptable() const { return !args_to_adapt_.empty(); }

    /// Return proposer scale.
    double scale() const { return is_adaptable() ? std::fabs((pDistr_->arg)(args_to_adapt_[0])) : hsp::na::DOUBLE; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Proposal : setter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    virtual void user_reset() {}
    virtual void associate_location(const Attribute &) {}

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Proposal : adaptation  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double get_accept_rate() const { return accept_rate_; }
    void update_accept_rate(bool is_accepted)
    {
        accept_rate_ = static_cast<double>(accept_rate_ * n_proposed_ + is_accepted) / (1. + n_proposed_);
        ++n_proposed_;
    }
    void reset_accept_rate()
    {
        accept_rate_ = 0.;
        n_proposed_ = 0;
    }
    /// Increase spread if acceptation rate is high and viceversa.
    void adapt(double ratio);

protected:
    std::string name_ = hsp::na::STRING;   /// name of the proposal
    std::shared_ptr<Distribution> pDistr_; /// distribution associated to proposal

    double user_scale_ = 1.;                         /// scale of the Metropolis proposal
    double accept_rate_ = 0.;                        /// acceptance rate
    unsigned long n_proposed_ = 0;                   /// n. of proposed moves
    std::vector<unsigned short> args_to_adapt_ = {}; /// which arguments should be adapted

    //===== propose : functional methods =====//
    double do_propose_sum_to_previous_(rnd_gen &prng, double xo) { return xo + pDistr_->rand(prng); }
    double do_propose_times_previous_(rnd_gen &prng, double xo) { return xo * pDistr_->rand(prng); }
    double do_propose_sample_rand_(rnd_gen &prng, double) { return pDistr_->rand(prng); }

    //===== logratio : functional methods =====//
    double do_logratio_pdf_(double xold, double xnew) const { return pDistr_->logratio(xold, xnew); }
    double do_logratio_symm_(double, double) const { return 0.; }
    double do_logratio_lognorm_symm_(double xold, double xnew) const { return std::log(xold / xnew); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Variable : class *** ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Variable : public Attribute
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Variable : constructor  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Variable(const std::string &name = hsp::na::STRING,
             const std::string &err_type = "variable",
             double value = hsp::na::DOUBLE,
             bool is_const = true,
             const hsp::parse::Definition &expr = hsp::na::STRING,
             bool is_prim = false)
        : pDistr_(std::make_unique<Distribution>()),
          pProp_(std::make_unique<Proposal>())
    {
        this->set_name(name, err_type);
        init(value, is_const, expr, is_prim);
    }
    ~Variable() = default;

    /// Sugar to call variable as a function if applicable, otherwise raise error.
    double operator()(double x) const { return (this->*evaluate_)(x); }

    /// Helper to initialize variable.
    void init(double value, bool is_const,
              const hsp::parse::Definition &parser = hsp::na::STRING, bool is_prim = false);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Variable : setter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Recover buffered value.
    void recover_value() { value_ = buffer_; }

    /// Restore user defined variable initialization.
    void user_reset(rnd_gen &prng);

    /// Assign an attribute to function.
    void init_func_attr(const std::string &key, const Attribute &attr) { pFunc_->init_attr(key, attr); }
    /// Initialize function reshaping. NOTE: to be called necessarily.
    void init_func_reshaping(bool discr);

    /// Initialize pointer to distribution.
    void init_distr(const hsp::parse::Definition &expr);
    /// Assign an attribute to distribution.
    void init_distr_attr(const std::string &key, const Attribute &attr);
    /// Finalize reshaping of distribution. NOTE: to be called necessarily.
    void init_distr_reshaping(bool discr);

    /// Random initialization of variable.
    void init_rand(rnd_gen &prng) { set(do_propose_sample_(prng)); }

    /// Initialize pointer to proposer.
    void init_proposer(const std::string &proposal = hsp::na::STRING,
                       double proposal_scale = hsp::na::DOUBLE);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Variable : getter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Get pointer to distribution.
    const std::shared_ptr<Distribution> get_distribution() const { return pDistr_; }

    bool is_constant() const override;
    bool is_discrete() const override final;
    bool is_model() const { return !pDistr_->is_empty(); }
    bool is_adaptable() const { return pProp_->is_adaptable(); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Variable : Metropolis-Hastings  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    void propose(rnd_gen &prng) { set((this->*propose_)(prng)); }
    bool is_valid_move() const { return (this->*is_valid_move_)(); }

    double get_distr_logpdf() const { return pDistr_->logpdf(value_); }
    double get_distr_logratio() const { return pDistr_->logratio(buffer_, value_); }
    double get_propose_logratio() const { return (this->*get_propose_logratio_)(); }
    double get_move_logratio() const { return (this->*get_move_logratio_)(); }

    //----- prop_ptr : adaptive methods -----//
    double get_accept_rate() const { return pProp_->get_accept_rate(); }
    void update_accept_rate(bool is_accepted) { pProp_->update_accept_rate(is_accepted); }
    void reset_accept_rate() { pProp_->reset_accept_rate(); }
    void adapt_proposal(double ratio) { pProp_->adapt(ratio); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Variable : show  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Get a description about the functional status.
    std::string str_summary() const override final;

    std::string distr_name() const { return pDistr_->name(); }
    std::string str_range() const { return pDistr_->str_range(); }
    std::string func_name() const { return pFunc_->name(); }
    std::string prop_name() const { return pProp_->name(); }

protected:
    double buffer_ = hsp::na::DOUBLE;     /// buffered value
    double user_value_ = hsp::na::DOUBLE; /// initial user defined value
    bool is_constant_ = true;             /// user flag for wheter constant or not

    // NOTE: Household must be copyable for simulations -> Variable pointers must be shared
    std::shared_ptr<Distribution> pDistr_; /// prior\model distribution of the variable
    std::shared_ptr<Proposal> pProp_;      /// Monte Carlo proposer

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Variable : getter [private]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double (Variable::*evaluate_)(double) const = &Variable::do_eval_error_;
    double do_eval_error_(double) const { throw std::runtime_error(
        "Called undefined evaluate operator for non-function " + err_name() + "."); }
    double do_eval_func_(double x) const { return pFunc_->evaluate(x); }

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

    //----- Propose : functional methods -----//
    double (Variable::*propose_)(rnd_gen &) = &Variable::do_propose_error_;
    double do_propose_error_(rnd_gen &) { throw std::invalid_argument(
        "Called to unknown proposer for " + err_name() + "."); }
    double do_propose_default_(rnd_gen &prng) { return pProp_->propose(prng, value_); }
    double do_propose_sample_(rnd_gen &prng) { return pDistr_->rand(prng); }

    //----- Valid move : functional methods -----//
    bool (Variable::*is_valid_move_)() const = &Variable::do_is_valid_move_default_;
    bool do_is_valid_move_default_() const { return pDistr_->is_valid(value_); }
    bool do_is_valid_move_sample_() const { return true; }

    /// Propose logratio - functional methods
    double (Variable::*get_propose_logratio_)() const = &Variable::do_get_propose_logratio_default_;
    double do_get_propose_logratio_default_() const { return pProp_->logratio(buffer_, value_); }
    double do_get_propose_logratio_sample_() const { return get_distr_logratio(); }

    //----- Move logratio : functional methods -----//
    double (Variable::*get_move_logratio_)() const = &Variable::do_get_move_logratio_default_;
    double do_get_move_logratio_default_() const { return get_distr_logratio() - (this->*get_propose_logratio_)(); }
    double do_move_logratio_sample_() const { return 0.; }
};

#endif // HOUSEPI__VARIABLE__H
