/*
 * variable.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "variable.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base_Proposer : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Base_Proposer::Base_Proposer(std::shared_ptr<Base_Distribution> distr,
                             double rule_start_value,
                             std::vector<int> args_to_adapt) : distr_ptr(distr)
{
    _args_to_adapt_ = args_to_adapt;

    if (rule_start_value == INF_DOUBLE)
    {
        _is_start_from_zero_ = false;
    }
    else
    {
        _is_start_from_zero_ = true;
        set_rule_start_value(rule_start_value);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - getter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Proposer::get_name() const { return distr_ptr->get_name(); }
double Base_Proposer::get_arg(int argID) const { return (distr_ptr->get_arg)(argID); }
double Base_Proposer::get_from_value() const { return from_value; }
bool Base_Proposer::is_rule_start_from_zero() const { return _is_start_from_zero_; }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Base_Proposer::set_rule_start_value(double x) { from_value = x; }
double Base_Proposer::rand(sys_rnd_gen &prng) { return (distr_ptr->rand)(prng); }
double Base_Proposer::propose(sys_rnd_gen &prng, double start) { return (this->*_func_propose_)(prng, start); }
double Base_Proposer::logratio(double xold, double xnew) const { return (this->*_logratio_method_)(xold, xnew); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Proposer - adaptation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Proposer::get_accept_rate() const
{
    return accept_rate;
}

/* Update the parameter accept rate. */
void Base_Proposer::update_accept_rate(bool is_accepted)
{
    accept_rate = (double)(accept_rate * n_proposed_moves + is_accepted) / (n_proposed_moves + 1);
    n_proposed_moves += 1;
}

/* Reset acceptation rate. */
void Base_Proposer::reset_accept_rate()
{
    accept_rate = 0.0;
    n_proposed_moves = 0;
};

/* Increase spread if acceptation rate is high and viceversa. */
void Base_Proposer::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);

    for (const int &argID : _args_to_adapt_)
    {
        (distr_ptr->set_arg)(argID, ((distr_ptr->get_arg)(argID)) * factor);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - show  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Base_Proposer::show_arg_names() const { return distr_ptr->show_arg_names(); }
const std::string Base_Proposer::show_arg_values() const { return distr_ptr->show_arg_values(); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Base_Proposer
{
public:
    //=====Constructors=====//
    Uniform_Proposer(
        dictargs rule_args = {{"min", (double)-0.5}, {"max", (double)0.5}},
        double rule_from = INF_DOUBLE) : Base_Proposer(distribution_switchboard("unif", rule_args), rule_from, {"min", "max"})
    {
        _func_propose_ = &Uniform_Proposer::_propose_sum_to_previous_;
        _logratio_method_ = &Uniform_Proposer::_logratio_symm_;
    }
    ~Uniform_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Normal_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Base_Proposer
{
public:
    //=====Constructor=====//
    Normal_Proposer(
        dictargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
        double rule_from = INF_DOUBLE) : Base_Proposer(distribution_switchboard("norm", rule_args), rule_from, {1})
    {
        _func_propose_ = &Normal_Proposer::_propose_sum_to_previous_;
        _logratio_method_ = is_rule_start_from_zero() ? &Normal_Proposer::_logratio_sum_to_zero_ : &Normal_Proposer::_logratio_symm_;
    }
    ~Normal_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Base_Proposer
{
public:
    //=====Constructor=====//
    Lognorm_Proposer(
        dictargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
        double rule_from = INF_DOUBLE) : Base_Proposer(distribution_switchboard("lnorm", rule_args), rule_from, {1})
    {
        _func_propose_ = &Lognorm_Proposer::_propose_times_previous_;
        _logratio_method_ = is_rule_start_from_zero() ? &Lognorm_Proposer::_logratio_times_zero_ : &Lognorm_Proposer::_logratio_lognorm_symm_;
    }
    ~Lognorm_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Gibbs_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gibbs_Proposer : public Base_Proposer
{
public:
    //=====Constructor=====//
    Gibbs_Proposer(std::shared_ptr<Base_Distribution> distr) : Base_Proposer(distr)
    {
        _func_propose_ = &Gibbs_Proposer::_propose_draw_rand_;
        _logratio_method_ = &Gibbs_Proposer::_logratio_pdf_;
    }
    ~Gibbs_Proposer() {}
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Switchboard to select the custom metropolis proposer rule. */
std::shared_ptr<Base_Proposer> proposer_switchboard(std::string rule_name, dictargs rule_args, double rule_from)
{
    std::shared_ptr<Base_Proposer> output = std::make_shared<Base_Proposer>(std::make_shared<Base_Distribution>());
    const std::vector<std::string> options = {"gibbs", "unif", "norm", "lnorm"};

    if (rule_name == "unif")
    {
        output = std::make_shared<Uniform_Proposer>(rule_args, rule_from);
    }
    else if (rule_name == "norm")
    {
        output = std::make_shared<Normal_Proposer>(rule_args, rule_from);
    }
    else if (rule_name == "lnorm")
    {
        output = std::make_shared<Lognorm_Proposer>(rule_args, rule_from);
    }
    else if (!tools::is_na(rule_name))
    {
        std::string err_message = "Unknown proposer rule '" + rule_name + "'.";
        err_message += tools::find_closest_option(rule_name, options);
        throw std::invalid_argument(err_message);
    }

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base_Variable : methods  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Base_Variable::assign_proposer(const std::string &rule_name,
                                    dictargs rule_args,
                                    double rule_from)
{
    // check that rule args have no mismatch
    if (rule_args.find(IO_MISMATCH) != rule_args.end())
    {
        tools::error("Mismatching definitions in rule arguments of '" + name + "'.");
    }

    if (tools::is_na(rule_name))
    {
        // immutable variable if no rule is defined
        infer = false;
        _func_propose_ = &Base_Variable::_propose_error_;
    }
    else
    {
        infer = true;

        // define proposer rule with functional programming
        if (rule_name == "gibbs")
        {
            prop_ptr = std::make_shared<Gibbs_Proposer>(distr_ptr);
            _func_propose_ = &Base_Variable::_propose_gibbs_sample_;
            _func_move_logratio_ = &Base_Variable::_move_logratio_gibbs_;
        }
        else
        {
            if ((rule_from != INF_DOUBLE) && (!distr_ptr->valid(rule_from)))
            {
                std::string min_val = distr_ptr->show_range_min();
                std::string max_val = distr_ptr->show_range_max();
                std::ostringstream oss;
                oss << "Initial value " << rule_from << " outside of value range"
                    << " {min:" << (tools::is_na(min_val) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(min_val))
                    << ", max:" << (tools::is_na(max_val) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(max_val))
                    << "}."
                    << " If intentional, ignore this message.";
                tools::warning(oss.str());
            }

            // assigne the pointer to the proposer
            prop_ptr = proposer_switchboard(rule_name, rule_args, rule_from);
            _func_propose_ = is_from_zero() ? &Base_Variable::_propose_from_zero_ : &Base_Variable::_propose_from_previous_;
        }
    }
}
