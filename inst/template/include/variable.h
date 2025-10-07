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
#include "parser.h"
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
             const std::vector<unsigned short> &args_to_adapt = {},
             bool is_discrete = false);
    Proposer(const std::string &distr_name,
             double rule_scale = HSP_NA_DOUBLE,
             const std::vector<unsigned short> &args_to_adapt = {},
             bool is_discrete = false);
    virtual ~Proposer() {};

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Proposer : getter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string name() const { return this->name_; }
    double get_arg(unsigned short argID) const { return pDistr_->arg(argID); }
    virtual double propose(rnd_gen &, double) { throw std::invalid_argument(
        "Called method of an undefined proposer."); }
    virtual double logratio(double, double) const { return 0.; }

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Proposer : setter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    virtual void user_reset() {}

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Proposer : adaptation  /////////////////////////////////
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

    double user_rule_scale_ = 1.; /// scale of the Metropolis rule

    double accept_rate_ = 0.;      /// acceptance rate
    unsigned long n_proposed_ = 0; /// n. of proposed moves

    unsigned short n_args_to_adapt_ = 0;           /// how many arguments to adapt
    unsigned short *args_to_adapt_list_ = nullptr; /// which arguments should be adapted

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

    Variable(const std::string &name = HSP_NA_STRING,
             const std::string &err_type = "variable",
             double value = HSP_NA_DOUBLE,
             bool is_const = true,
             bool is_discr = false,
             const DefinitionParser &expr = HSP_NA_STRING,
             bool is_prim = false)
        : pDistr_(std::make_shared<Distribution>()),
          pProp_(std::make_shared<Proposer>())
    {
        this->name(name, err_type);
        init(value, is_const, is_discr, expr, is_prim);
    }
    ~Variable() = default;

    /// Sugar to call variable as a function if applicable, otherwise raise error.
    double operator()(double x) const { return (this->*evaluate_)(x); }

    /// Helper to initialize variable.
    void init(double value,
              bool is_const, bool is_discr,
              const DefinitionParser &parser, bool is_prim);

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
    void init_func_reshaping();

    /// Initialize pointer to distribution.
    void init_distr(const DefinitionParser &expr);
    /// Assign an attribute to distribution.
    void init_distr_attr(const std::string &key, const Attribute &attr);
    /// Finalize reshaping of distribution. NOTE: to be called necessarily.
    void init_distr_reshaping();

    /// Initialize location parameter.
    void init_location_arg(const Attribute &x);

    /// Random initialization of variable.
    void init_rand(rnd_gen &prng);

    /// Initialize pointer to proposer.
    void init_proposer(const std::string &rule_name = HSP_NA_STRING,
                       double rule_scale = HSP_NA_DOUBLE);

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  Variable : getter  ///////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    bool is_constant() const override ;
    bool is_discrete() const override final;
    bool is_model() const { return !pDistr_->is_empty(); }
    /// Get a description about the functional status.
    std::string str_summary() const override final;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Variable : Metropolis-Hastings  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    void propose(rnd_gen &prng) { set((this->*propose_)(prng)); }
    bool is_valid_move() const { return (this->*is_valid_move_)(); }

    double get_distr_logpdf() const { return (this->*distr_logpdf_)(); }
    double get_distr_logratio() const { return (this->*distr_logratio_)(); }
    double get_propose_logratio() const { return (this->*get_propose_logratio_)(); }
    double get_move_logratio() const { return (this->*get_move_logratio_)(); }

    //----- prop_ptr : adaptive methods -----//
    double get_accept_rate() const { return pProp_->get_accept_rate(); }
    void update_accept_rate(bool is_accepted) { pProp_->update_accept_rate(is_accepted); }
    void reset_accept_rate() { pProp_->reset_accept_rate(); }
    void adapt_rule(double ratio) { pProp_->adapt(ratio); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  Variable : show  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::string distr_name() const { return pDistr_->name(); }
    std::string distr_str_arg_names() const { return pDistr_->str_arg_names(); }
    std::string distr_str_arg_values() const { return pDistr_->str_arg_values(); }

    std::string str_min() const { return pDistr_->str_min(); }
    std::string str_max() const { return pDistr_->str_max(); }
    std::string str_range() const { return pDistr_->str_range(); }

    std::string func_name() const { return pFunc_->name(); }

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
    double buffer_ = HSP_NA_DOUBLE;     /// buffered value
    double user_value_ = HSP_NA_DOUBLE; /// initial user defined value
    Attribute loc_;                     /// location parameter
    bool is_constant_ = true;           /// user flag for wheter constant or not
    bool is_discrete_ = false;          /// user flag for wheter discrete or not

    std::shared_ptr<Distribution> pDistr_; /// prior\model distribution of the variable
    std::shared_ptr<Proposer> pProp_;      /// Monte Carlo proposer

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Variable : getter [private]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double (Variable::*evaluate_)(double) const = &Variable::do_eval_error_;
    double do_eval_error_(double) const { throw std::runtime_error(
        "Wrong call to evaluate operator for non-function " + err_name() + "."); }
    double do_eval_func_(double x) const { return pFunc_->evaluate(x); }

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  Variable : setter [private]  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Helper with proposer initialization of direct sample from prior/model.
    void helper_init_sample_proposer_();

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
    double do_propose_with_loc_(rnd_gen &prng) { return loc_.get() + pProp_->propose(prng, value_ - loc_.get()); }
    double do_propose_sample_(rnd_gen &prng) { return pDistr_->rand(prng); }
    double do_propose_sample_with_loc_(rnd_gen &prng) { return loc_.get() + do_propose_sample_(prng); }

    //----- Valid move : functional methods -----//
    bool (Variable::*is_valid_move_)() const = &Variable::do_is_valid_move_default_;
    bool do_is_valid_move_default_() const { return pDistr_->is_valid(value_); }
    bool do_is_valid_move_with_loc_() const { return pDistr_->is_valid(value_ - loc_.get()); }
    bool do_is_valid_move_sample_() const { return true; }

    ///  Distribution logpdf - function methods
    double (Variable::*distr_logpdf_)() const = &Variable::do_distr_logpdf_default_;
    double do_distr_logpdf_default_() const { return pDistr_->logpdf(value_); }
    double do_distr_logpdf_with_loc_() const { return pDistr_->logpdf(value_ - loc_.get()); }

    ///  Distribution logpdf - function methods
    double (Variable::*distr_logratio_)() const = &Variable::do_distr_logratio_default_;
    double do_distr_logratio_default_() const { return pDistr_->logratio(buffer_, value_); }
    double do_distr_logratio_with_loc_() const
    {
        const double loc_val = loc_.get();
        return pDistr_->logratio(buffer_ - loc_val, value_ - loc_val);
    }

    /// Propose logratio - functional methods
    double (Variable::*get_propose_logratio_)() const = &Variable::do_get_propose_logratio_default_;
    double do_get_propose_logratio_default_() const { return pProp_->logratio(buffer_, value_); }
    double do_get_propose_logratio_with_loc_() const
    {
        const double loc_val = loc_.get();
        return pProp_->logratio(buffer_ - loc_val, value_ - loc_val);
    }
    double do_get_propose_logratio_sample_() const { return do_distr_logratio_default_(); }
    double do_get_propose_logratio_sample_with_loc_() const { return do_distr_logratio_with_loc_(); }

    //----- Move logratio : functional methods -----//
    double (Variable::*get_move_logratio_)() const = &Variable::do_get_move_logratio_default_;
    double do_get_move_logratio_default_() const
    {
        return pDistr_->logratio(buffer_, value_) - pProp_->logratio(buffer_, value_);
    }
    double do_get_move_logratio_with_loc_() const
    {
        const double loc_val = loc_.get();
        const double xo = buffer_ - loc_val, xn = value_ - loc_val;
        return pDistr_->logratio(xo, xn) - pProp_->logratio(xo, xn);
    }
    double do_move_logratio_sample_() const { return 0.; }
};

#endif // HOUSEPI__VARIABLE__H
