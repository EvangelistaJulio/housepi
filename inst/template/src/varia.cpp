/*
 * varia.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "function.h"
#include "varia.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : methods  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double _Base_Distribution_::pdf(double x) const
{
    return (this->*_func_trunc_)(_pdf_(x));
}
double _Base_Distribution_::cdf(double x) const
{
    return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(x)));
}
double _Base_Distribution_::logpdf(double x) const
{
    return (this->*_func_logtrunc_)(_logpdf_(x));
}

double _Base_Distribution_::logratio(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.0;
}
double _Base_Distribution_::rand(sys_rnd_gen &prng)
{
    UNUSED(prng);
    return _undefined_call_error_("rand");
}

double _Base_Distribution_::_pdf_(double x) const { return _undefined_call_error_("pdf", x); }
double _Base_Distribution_::_cdf_(double x) const { return _undefined_call_error_("cdf", x); }
double _Base_Distribution_::_logpdf_(double x) const { return _undefined_call_error_("logpdf", x); }

/* Internal error for calls to undefined distribution. */
double _Base_Distribution_::_undefined_call_error_(const std::string &txt, double x, double y) const
{
    UNUSED(x);
    UNUSED(y);
    tools::error("Attempt to call undefined method '" + txt + "' of distribution '" + get_name() + "'.");
    return std::numeric_limits<double>::quiet_NaN();
}

/* Truncate distribution at declaration.  */
void _Base_Distribution_::set_truncation(const dictargs &range)
{
    int trunc_cases = set_range(range);
    if (trunc_cases == 0) // no truncation
    {
        _func_trunc_ = &_Base_Distribution_::_trunc_dummy_;
        _func_logtrunc_ = &_Base_Distribution_::_logtrunc_dummy_;
        _func_shift_ = &_Base_Distribution_::_shift_dummy_;
    }
    else if (trunc_cases == 1) // lower range truncation
    {
        _func_trunc_ = &_Base_Distribution_::_trunc_min_;
        _func_logtrunc_ = &_Base_Distribution_::_logtrunc_min_;
        _func_shift_ = &_Base_Distribution_::_shift_min_;
    }
    else if (trunc_cases == 2) // upper range truncation
    {
        _func_trunc_ = &_Base_Distribution_::_trunc_max_;
        _func_logtrunc_ = &_Base_Distribution_::_logtrunc_max_;
        _func_shift_ = &_Base_Distribution_::_shift_dummy_;
    }
    else // both sides truncation
    {
        _func_trunc_ = &_Base_Distribution_::_trunc_both_sides_;
        _func_logtrunc_ = &_Base_Distribution_::_logtrunc_both_sides_;
        _func_shift_ = &_Base_Distribution_::_shift_dummy_;
    }
}

/* Error for. */
void _Base_Distribution_::check_rand_attempts(int &num_attempts) const
{
    if (num_attempts > MAX_NUMBER_ATTEMPTS)
    {
        std::ostringstream oss;
        oss << "Exceeded attempts (N > " << MAX_NUMBER_ATTEMPTS << ")"
            << " in random draw of " << get_name() << "(" << show_arg_values() << ")"
            << " with truncation {" << get_range_min() << ", " << get_range_max() << "}.";
        tools::error(oss.str());
    }
    ++num_attempts;
}

//----- Functional methods for the truncation -----//
double _Base_Distribution_::_trunc_dummy_(double x) const { return x; }
double _Base_Distribution_::_trunc_min_(double x) const { return x / (1.0 - _cdf_(get_range_min())); }
double _Base_Distribution_::_trunc_max_(double x) const { return x / (_cdf_(get_range_max())); }
double _Base_Distribution_::_trunc_both_sides_(double x) const { return x / (_cdf_(get_range_max()) - _cdf_(get_range_max())); }

//----- Functional methods for the (log) truncation -----//
double _Base_Distribution_::_logtrunc_dummy_(double x) const { return x; }
double _Base_Distribution_::_logtrunc_min_(double x) const { return x - log(1.0 - _cdf_(get_range_min())); }
double _Base_Distribution_::_logtrunc_max_(double x) const { return x - log(_cdf_(get_range_max())); }
double _Base_Distribution_::_logtrunc_both_sides_(double x) const { return x - log(_cdf_(get_range_max()) - _cdf_(get_range_min())); }

//----- Functional methods to shift cumulative -----//
double _Base_Distribution_::_shift_dummy_(double x) const { return x; }
double _Base_Distribution_::_shift_min_(double x) const { return x - _cdf_(get_range_min()); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  _Base_Proposer_ : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

_Base_Proposer_::_Base_Proposer_(std::shared_ptr<_Base_Distribution_> distr,
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
//////////////////////////////////  _Base_Proposer_ : methods  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Increase spread if acceptation rate is high and viceversa. */
void _Base_Proposer_::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);

    for (const int &argID : _args_to_adapt_)
    {
        (distr_ptr->set_arg)(argID, ((distr_ptr->get_arg)(argID)) * factor);
    }
}

/* Reset acceptation rate. */
void _Base_Proposer_::reset_accept_rate()
{
    accept_rate = 0.0;
    num_proposed_moves = 0;
};

/* Update the parameter accept rate. */
void _Base_Proposer_::update_accept_rate(bool is_accepted)
{
    accept_rate = (double)(accept_rate * num_proposed_moves + is_accepted) / (num_proposed_moves + 1);
    num_proposed_moves += 1;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public _Base_Proposer_
{
public:
    //=====Constructors=====//
    Uniform_Proposer(
        dictargs rule_args = {{"min", (double)-0.5}, {"max", (double)0.5}},
        double rule_from = INF_DOUBLE) : _Base_Proposer_(distribution_ptr("unif", rule_args), rule_from, {"min", "max"})
    {
        _func_propose_ = &Uniform_Proposer::_propose_sum_to_previous_;
        _logratio_method_ = &Uniform_Proposer::_logratio_symm_;
    }
    ~Uniform_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Normal_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public _Base_Proposer_
{
public:
    //=====Constructor=====//
    Normal_Proposer(
        dictargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
        double rule_from = INF_DOUBLE) : _Base_Proposer_(distribution_ptr("norm", rule_args), rule_from, {1})
    {
        _func_propose_ = &Normal_Proposer::_propose_sum_to_previous_;
        _logratio_method_ = is_rule_start_from_zero() ? &Normal_Proposer::_logratio_sum_to_zero_ : &Normal_Proposer::_logratio_symm_;
    }
    ~Normal_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public _Base_Proposer_
{
public:
    //=====Constructor=====//
    Lognorm_Proposer(
        dictargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
        double rule_from = INF_DOUBLE) : _Base_Proposer_(distribution_ptr("lnorm", rule_args), rule_from, {1})
    {
        _func_propose_ = &Lognorm_Proposer::_propose_times_previous_;
        _logratio_method_ = is_rule_start_from_zero() ? &Lognorm_Proposer::_logratio_times_zero_ : &Lognorm_Proposer::_logratio_lognorm_symm_;
    }
    ~Lognorm_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Gibbs_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gibbs_Proposer : public _Base_Proposer_
{
public:
    //=====Constructor=====//
    Gibbs_Proposer(std::shared_ptr<_Base_Distribution_> distr) : _Base_Proposer_(distr)
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
std::shared_ptr<_Base_Proposer_> proposer_switchboard(std::string rule_name, dictargs rule_args, double rule_from)
{
    std::shared_ptr<_Base_Proposer_> output = std::make_shared<_Base_Proposer_>(std::make_shared<_Base_Distribution_>());

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
        throw std::invalid_argument(
            "Unknown proposer rule '" + rule_name + "' for Metropolis-Hastings.");
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
