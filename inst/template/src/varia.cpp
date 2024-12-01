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
////////////////////////////////  Base Distribution : virtual  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_pdf_(double x) const { return _undefined_call_error_("pdf", x); }
double Base_Distribution::_cdf_(double x) const { return _undefined_call_error_("cdf", x); }
double Base_Distribution::_logpdf_(double x) const { return _undefined_call_error_("logpdf", x); }
double Base_Distribution::_logratio_(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.0;
}
double Base_Distribution::_rand_(sys_rnd_gen &prng)
{
    UNUSED(prng);
    return _undefined_call_error_("rand");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : deaults  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_default_pdf_(double x) const { return (this->*_func_trunc_)(_pdf_(x)); }
double Base_Distribution::_default_cdf_(double x) const { return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(x))); }
double Base_Distribution::_default_logpdf_(double x) const { return (this->*_func_logtrunc_)(_logpdf_(x)); }
double Base_Distribution::_default_logratio_(double xold, double xnew) const { return _logratio_(xold, xnew); }
double Base_Distribution::_default_rand_(sys_rnd_gen &prng) { return _rand_(prng); }

/* Internal error for calls to undefined distribution. */
double Base_Distribution::_undefined_call_error_(const std::string &txt, double x, double y) const
{
    UNUSED(x);
    UNUSED(y);
    tools::error("Call to undefined method '" + txt + "' of distribution '" + get_name() + "'.");
    return std::numeric_limits<double>::quiet_NaN();
}

/* Check accept/reject within range is not exceeding maximum number of attempts. */
void Base_Distribution::_check_rand_attempts_(int &n_attempts) const
{
    if (n_attempts > MAX_NUMBER_ATTEMPTS)
    {
        std::ostringstream oss;
        oss << "Exceeded attempts (N > " << MAX_NUMBER_ATTEMPTS << ")"
            << " in random draw of " << get_name() << "(" << show_arg_values() << ")"
            << " with truncation {" << get_range_min() << ", " << get_range_max() << "}.";
        tools::error(oss.str());
    }
    ++n_attempts;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized center /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_discr_center_pdf_(double x) const
{
    double d = std::round(x);
    return (this->*_func_trunc_)(_cdf_(d + 0.5) - _cdf_(d - 0.5));
}

double Base_Distribution::_discr_center_cdf_(double x) const
{
    double d = std::round(x);
    return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(d + 0.5)));
}

double Base_Distribution::_discr_center_logpdf_(double x) const
{
    return log(_discr_center_pdf_(x));
}

double Base_Distribution::_discr_center_logratio_(double xold, double xnew) const
{
    double dold = std::round(xold), dnew = std::round(xnew);
    return log((_cdf_(dnew + 0.5) - _cdf_(dnew - 0.5)) / (_cdf_(dold + 0.5) - _cdf_(dold - 0.5)));
}

double Base_Distribution::_discr_center_rand_(sys_rnd_gen &prng)
{
    return std::round(_rand_(prng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : discretized left /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_discr_left_pdf_(double x) const
{
    double d = std::floor(x);
    return (this->*_func_trunc_)(_cdf_(d + 1.0) - _cdf_(d));
}

double Base_Distribution::_discr_left_cdf_(double x) const
{
    double d = std::floor(x);
    return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(d + 1.0)));
}

double Base_Distribution::_discr_left_logpdf_(double x) const
{
    return log(_discr_left_pdf_(x));
}

double Base_Distribution::_discr_left_logratio_(double xold, double xnew) const
{
    double dold = std::floor(xold), dnew = std::floor(xnew);
    return log((_cdf_(dnew + 1.0) - _cdf_(dnew)) / (_cdf_(dold + 1.0) - _cdf_(dold)));
}

double Base_Distribution::_discr_left_rand_(sys_rnd_gen &prng)
{
    return std::floor(_rand_(prng));
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Base Distribution : truncation (functional)  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_trunc_dummy_(double x) const { return x; }
double Base_Distribution::_trunc_min_(double x) const { return x / (1.0 - _cdf_(get_range_min())); }
double Base_Distribution::_trunc_max_(double x) const { return x / (_cdf_(get_range_max())); }
double Base_Distribution::_trunc_both_sides_(double x) const { return x / (_cdf_(get_range_max()) - _cdf_(get_range_max())); }

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Base Distribution : logtruncation (functional)  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_logtrunc_dummy_(double x) const { return x; }
double Base_Distribution::_logtrunc_min_(double x) const { return x - log(1.0 - _cdf_(get_range_min())); }
double Base_Distribution::_logtrunc_max_(double x) const { return x - log(_cdf_(get_range_max())); }
double Base_Distribution::_logtrunc_both_sides_(double x) const { return x - log(_cdf_(get_range_max()) - _cdf_(get_range_min())); }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  Base Distribution : shift cumulative (functional)  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_shift_dummy_(double x) const { return x; }
double Base_Distribution::_shift_min_(double x) const { return x - _cdf_(get_range_min()); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : assign reshaping  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Truncate and/or discretize distribution at declaration.  */
void Base_Distribution::_assign_reshaping_(const dictargs &range, const std::string &discr)
{
    dictargs new_range = range;

    //===== DISCRETIZATION =====//
    // asjust range in case of discretization
    if (discr != UNDEFINED)
    {
        if (discr == "center")
        { // centered discretization
            _func_pdf_ = &Base_Distribution::_discr_center_pdf_;
            _func_cdf_ = &Base_Distribution::_discr_center_cdf_;
            _func_logpdf_ = &Base_Distribution::_discr_center_logpdf_;
            _func_logratio_ = &Base_Distribution::_discr_center_logratio_;
            _func_rand_ = &Base_Distribution::_discr_center_rand_;

            if (range.find("min") != range.end())
            {
                new_range["min"] = std::floor(range.at("min")) + 0.5;
            }
            if (range.find("max") != range.end())
            {
                new_range["max"] = std::ceil(range.at("max")) - 0.5;
            }
        }
        else if (discr == "left")
        { // left adjusted discretization
            _func_pdf_ = &Base_Distribution::_discr_left_pdf_;
            _func_cdf_ = &Base_Distribution::_discr_left_cdf_;
            _func_logpdf_ = &Base_Distribution::_discr_left_logpdf_;
            _func_logratio_ = &Base_Distribution::_discr_left_logratio_;
            _func_rand_ = &Base_Distribution::_discr_left_rand_;
            
            if (range.find("min") != range.end())
            {
                new_range["min"] = std::floor(range.at("min"));
            }
            if (range.find("max") != range.end())
            {
                new_range["max"] = std::ceil(range.at("max"));
            }
        }
        else
        {
            tools::error(
                "Unrecognized discretization method '" + discr + "'.");
        }
    }

    //===== TRUNCATION =====/
    int trunc_cases = set_range(new_range);

    //----- Fucntional methods -----//
    if (trunc_cases == 0) // no truncation
    {
        _func_trunc_ = &Base_Distribution::_trunc_dummy_;
        _func_logtrunc_ = &Base_Distribution::_logtrunc_dummy_;
        _func_shift_ = &Base_Distribution::_shift_dummy_;
    }
    else if (trunc_cases == 1) // lower range truncation
    {
        _func_trunc_ = &Base_Distribution::_trunc_min_;
        _func_logtrunc_ = &Base_Distribution::_logtrunc_min_;
        _func_shift_ = &Base_Distribution::_shift_min_;
    }
    else if (trunc_cases == 2) // upper range truncation
    {
        _func_trunc_ = &Base_Distribution::_trunc_max_;
        _func_logtrunc_ = &Base_Distribution::_logtrunc_max_;
        _func_shift_ = &Base_Distribution::_shift_dummy_;
    }
    else // both sides truncation
    {
        _func_trunc_ = &Base_Distribution::_trunc_both_sides_;
        _func_logtrunc_ = &Base_Distribution::_logtrunc_both_sides_;
        _func_shift_ = &Base_Distribution::_shift_dummy_;
    }
}

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
//////////////////////////////////  Base_Proposer : methods  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Increase spread if acceptation rate is high and viceversa. */
void Base_Proposer::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);

    for (const int &argID : _args_to_adapt_)
    {
        (distr_ptr->set_arg)(argID, ((distr_ptr->get_arg)(argID)) * factor);
    }
}

/* Reset acceptation rate. */
void Base_Proposer::reset_accept_rate()
{
    accept_rate = 0.0;
    n_proposed_moves = 0;
};

/* Update the parameter accept rate. */
void Base_Proposer::update_accept_rate(bool is_accepted)
{
    accept_rate = (double)(accept_rate * n_proposed_moves + is_accepted) / (n_proposed_moves + 1);
    n_proposed_moves += 1;
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

class Uniform_Proposer : public Base_Proposer
{
public:
    //=====Constructors=====//
    Uniform_Proposer(
        dictargs rule_args = {{"min", (double)-0.5}, {"max", (double)0.5}},
        double rule_from = INF_DOUBLE) : Base_Proposer(distribution_ptr("unif", rule_args), rule_from, {"min", "max"})
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
        double rule_from = INF_DOUBLE) : Base_Proposer(distribution_ptr("norm", rule_args), rule_from, {1})
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
        double rule_from = INF_DOUBLE) : Base_Proposer(distribution_ptr("lnorm", rule_args), rule_from, {1})
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
