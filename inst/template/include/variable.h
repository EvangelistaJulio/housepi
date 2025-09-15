/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * variable.h
 *
 */

#ifndef HOUSEPI__VARIABLE__H
#define HOUSEPI__VARIABLE__H

#include <memory>
#include <string>
#include <vector>

#include "tools.h"
#include "distribution.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// *** Proposer : class *** ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Base class for Metropolis proposers.
class Proposer
{
public:
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  Proposer : constructor  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    Proposer() {};
    Proposer(std::shared_ptr<Distribution>,
             double rule_scale = HSP_NA_DOUBLE,
             double rule_start_value = HSP_NA_DOUBLE,
             const std::vector<unsigned short> &args_to_adapt = {},
             bool is_discrete = false);
    Proposer(const std::string &distr_name,
             double rule_scale = HSP_NA_DOUBLE,
             double rule_start_value = HSP_NA_DOUBLE,
             const std::vector<unsigned short> &args_to_adapt = {},
             bool is_discrete = false);
    virtual ~Proposer() {};

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Proposer : getter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string name() const { return this->name_; }
    double get_arg(unsigned short argID) const { return (pDistr_->arg)(argID); }
    double get_from_value() const { return from_value_; }
    double propose(rnd_gen &prng, double start) { return (this->*propose_)(prng, start); }
    double logratio(double xold, double xnew) const { return (this->*logratio_)(xold, xnew); }
    double rand(rnd_gen &prng) { return (pDistr_->rand)(prng); }
    bool is_rule_start_from_zero() const { return is_start_from_zero_; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Proposer : setter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    virtual void user_reset() {}
    void set_rule_start_value(double x) { from_value_ = x; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Proposer : adaptation  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double get_accept_rate() const { return accept_rate_; }
    void update_accept_rate(bool is_accepted)
    {
        accept_rate_ = static_cast<double>(accept_rate_ * n_proposed_moves_ + is_accepted) / (1. + n_proposed_moves_);
        ++n_proposed_moves_;
    }
    void reset_accept_rate()
    {
        accept_rate_ = 0.;
        n_proposed_moves_ = 0;
    }
    void adapt(double ratio);

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Proposer : show  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string str_arg_names() const { return pDistr_->str_arg_names(); }
    std::string str_arg_values() const { return pDistr_->str_arg_values(); }
    std::string str_scale() const;

protected:
    std::string name_ = HSP_NA_STRING;
    std::shared_ptr<Distribution> pDistr_;

    double user_rule_scale_ = 1.;       /// scale of the Metropolis rule
    bool is_start_from_zero_ = false;   /// wheter to move from given value
    double from_value_ = HSP_NA_DOUBLE; /// given value where to start from

    double accept_rate_ = 0.;            /// acceptance rate
    unsigned long n_proposed_moves_ = 0; /// n. of proposed moves

    unsigned short n_args_to_adapt_ = 0;           /// how many arguments to adapt
    unsigned short *args_to_adapt_list_ = nullptr; /// which arguments should be adapted

    //===== propose : functional methods =====//
    double (Proposer::*propose_)(rnd_gen &, double) = &Proposer::do_propose_undefined_;
    double do_propose_undefined_(rnd_gen &, double) { throw std::invalid_argument(
        "Called method of an undefined proposer."); }
    double do_propose_sum_to_previous_(rnd_gen &prng, double start) { return start + rand(prng); }
    double do_propose_times_previous_(rnd_gen &prng, double start) { return start * rand(prng); }
    double do_propose_sample_rand_(rnd_gen &prng, double) { return rand(prng); }

    //===== logratio : functional methods =====//
    double (Proposer::*logratio_)(double, double) const = &Proposer::do_logratio_undefined_;
    double do_logratio_undefined_(double, double) const { return 0.; }
    double do_logratio_pdf_(double xold, double xnew) const { return pDistr_->logratio(xold, xnew); }
    double do_logratio_sum_to_zero_(double xold, double xnew) const { return do_logratio_pdf_(xold - from_value_, xnew - from_value_); }
    double do_logratio_times_zero_(double xold, double xnew) const { return do_logratio_pdf_(xold / from_value_, xnew / from_value_); }
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

    Variable(const std::string &name = HSP_NA_STRING,
             double value = HSP_NA_DOUBLE,
             bool is_const = true,
             bool is_discr = false,
             const std::string &func_name = HSP_NA_STRING,
             bool is_prim = false)
        : pDistr_(std::make_shared<Distribution>()),
          pProp_(std::make_shared<Proposer>())
    {
		err_type_ = "variable";
        init(name, value, is_const, is_discr, func_name, is_prim);
    }
    ~Variable() = default;

    /// Sugar to call variable as a function if applicable, otherwise raise error.
    double operator()(double x) const { return (this->*evaluate_)(x); }

    /// Helper to initialize variable.
    void init(const std::string &name, double value,
              bool is_const = true, bool is_discr = false,
              const std::string &func_name = HSP_NA_STRING,
              bool is_prim = false);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Variable : setter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Recover buffered value.
    void recover_value() { value_ = buffer_value_; }

    /// Restore user defined variable initialization.
    void user_reset()
    {
        value_ = user_value_;
        pProp_->user_reset();
    }

    /// Assign a scalar argument to function.
    void init_func_attr(const std::string &a, bool is_arg, double x) { pFunc_->init_attr(a, is_arg, x); }
    /// Assign a pointer argument to function.
    void init_func_attr(const std::string &a, bool is_arg, observ_dbl_ptr y) { pFunc_->init_attr(a, is_arg, y); }
    /// Assign a function attribute to function.
    void init_func_attr(const std::string &a, bool is_arg, std::shared_ptr<Function> f) { pFunc_->init_attr(a, is_arg, f); }
    /// Initialize function reshaping.
    void init_func_reshaping();

    /// Initialize pointer to distribution.
    void init_distr(const std::string &distr_name = HSP_NA_STRING);
    /// Assign a scalar argument to distribution.
    void init_distr_attr(const std::string &a, bool is_arg, double x) { pDistr_->init_attr(a, is_arg, x); }
    /// Assign a pointer argument to distribution.
    void init_distr_attr(const std::string &a, bool is_arg, observ_dbl_ptr y) { pDistr_->init_attr(a, is_arg, y); }
    /// Assign a function argument to distribution.
    void init_distr_attr(const std::string &a, bool is_arg, std::shared_ptr<Function> f) { pDistr_->init_attr(a, is_arg, f); }
    /// Finalize reshaping of distribution.
    void init_distr_reshaping();

    /// Initialize pointer to proposer.
    void init_proposer(const std::string &rule_name = HSP_NA_STRING,
                       double rule_scale = HSP_NA_DOUBLE,
                       double rule_from = HSP_NA_DOUBLE);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Variable : getter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    bool is_constant() const override final { return Attribute::is_constant() && is_constant_; }
    bool is_discrete() const { return pDistr_->is_discrete(); }
    bool is_from_zero() const { return (pProp_->is_rule_start_from_zero)(); }
    double get_from_value() const { return (pProp_->get_from_value)(); }
    bool is_chain_dependency(const Variable &other) const;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Variable : Metropolis-Hastings  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double get_move_logratio() const { return (this->*move_logratio_)(); }
    bool is_valid_move() const { return (this->*is_valid_move_)(); }

    //----- distr_ptr : prior/model methods -----//
    void rand(rnd_gen &prng) { set(pDistr_->rand(prng)); }

    double distr_logpdf() const { return pDistr_->logpdf(value_); }
    double distr_logratio() const { return pDistr_->logratio(buffer_value_, value_); }

    //----- prop_ptr : propose methods -----//
    void propose(rnd_gen &prng) { set((this->*propose_)(prng)); }
    double propose_logratio() const { return (pProp_->logratio)(buffer_value_, value_); }

    //----- prop_ptr : adaptive methods -----//
    double get_accept_rate() const { return pProp_->get_accept_rate(); }
    void update_accept_rate(bool is_accepted) { pProp_->update_accept_rate(is_accepted); }
    void reset_accept_rate() { pProp_->reset_accept_rate(); }
    void adapt_rule(double ratio) { (pProp_->adapt)(ratio); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Variable : show  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string distr_name() const { return pDistr_->name(); }
    std::string distr_str_arg_names() const { return pDistr_->str_arg_names(); }
    std::string distr_str_arg_values() const { return pDistr_->str_arg_values(); }
    std::string distr_str_summary() const { return pDistr_->str_summary(); }

    std::string str_min() const { return pDistr_->str_min(); }
    std::string str_max() const { return pDistr_->str_max(); }

    std::string func_name() const { return pFunc_->name(); }
    std::string func_str_summary() const { return pFunc_->str_summary(); }

    std::string prop_name() const { return pProp_->name(); }
    std::string prop_str_arg_names() const { return pProp_->str_arg_names(); }
    std::string prop_str_arg_values() const { return pProp_->str_arg_values(); }
    std::string prop_str_scale() const { return pProp_->str_scale(); }

#ifdef _HOUSEPI_RCPP
    Rcpp::NumericVector rcpp_export_distr_args() const { return pDistr_->rcpp_export_args(); }
    Rcpp::NumericVector rcpp_export_func_args() const { return pFunc_->rcpp_export_args(); }
    Rcpp::NumericVector rcpp_export_range() const { return pDistr_->rcpp_export_range(); }
#endif

protected:
    double buffer_value_ = HSP_NA_DOUBLE; /// buffered value
    double user_value_ = HSP_NA_DOUBLE;   /// initial user defined value
    bool is_discrete_ = false;            /// user flag for wheter discrete or not
    bool is_constant_ = false;            /// user flag for wheter constant or not

    std::shared_ptr<Distribution> pDistr_; /// prior\model distribution of the variable
    std::shared_ptr<Proposer> pProp_;      /// Monte Carlo proposer

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Variable : getter [private]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double (Variable::*evaluate_)(double) const = &Variable::do_eval_error_;
    double do_eval_error_(double) const { throw std::runtime_error(
        "Wrong call to evaluate operator for non-function " + err_name() + "."); }
    double do_eval_func_(double x) const { return pFunc_->evaluate(x); }

    double do_get_default_() const { return value_; }
    double do_get_func_eval_() const { return pFunc_->evaluate(); }

    observ_dbl_ptr do_get_ptr_default_() const { return &value_; }
    observ_dbl_ptr do_get_ptr_func_eval_() const { throw std::runtime_error(
        "Attempt to call pointer to function \"" + name() + "\"."); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Variable : setter [private]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Store the buffer and change value.
    void do_set_default_(double x) override final
    {
        buffer_value_ = value_;
        value_ = x;
    }
    /// Set error.
    void do_set_error_(double) override final { throw std::runtime_error(
        "Non scalar " + err_name() + " cannot be changed."); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  Variable : Metropolis-Hastings [private]  ///////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    //----- Valid move : functional methods -----//
    bool (Variable::*is_valid_move_)() const = &Variable::do_is_valid_move_default_;
    bool do_is_valid_move_default_() const { return pDistr_->is_valid(value_); }
    bool do_is_valid_move_sample_() const { return true; }

    //----- Move logratio : functional methods -----//
    double (Variable::*move_logratio_)() const = &Variable::do_move_logratio_default_;
    double do_move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
    double do_move_logratio_sample_() const { return 0.; }

    //----- Propose : functional methods -----//
    double (Variable::*propose_)(rnd_gen &) = &Variable::do_propose_error_;
    double do_propose_error_(rnd_gen &) { throw std::invalid_argument(
        "Called to unknown proposer for " + err_name() + "."); }
    double do_propose_sample_(rnd_gen &prng) { return (pProp_->propose)(prng, 0); }
    double do_propose_from_previous_(rnd_gen &prng) { return (pProp_->propose)(prng, value_); }
    double do_propose_from_zero_(rnd_gen &prng) { return (pProp_->propose)(prng, pProp_->get_from_value()); }
};

#endif // HOUSEPI__VARIABLE__H
