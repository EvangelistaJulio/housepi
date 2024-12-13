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
                             std::vector<ushort> args_to_adapt) : distr_ptr(distr)
{
    // by default proposal is from previous state
    _is_start_from_zero_ = false; 

    // proposal alwats starting by given value
    if (rule_start_value != INF_DOUBLE)
    {
        _is_start_from_zero_ = true;
        set_rule_start_value(rule_start_value);
    }

    // load which are the parameter to be adapted
    _args_to_adapt_ = args_to_adapt;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - getter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Proposer::get_name() const { return distr_ptr->get_name(); }
double Base_Proposer::get_arg(ushort argID) const { return (distr_ptr->get_arg)(argID); }
double Base_Proposer::get_from_value() const { return from_value; }
bool Base_Proposer::is_rule_start_from_zero() const { return _is_start_from_zero_; }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Base_Proposer::set_rule_start_value(double x) { from_value = x; }
double Base_Proposer::rand(rnd_gen &prng) { return (distr_ptr->rand)(prng); }
double Base_Proposer::propose(rnd_gen &prng, double start) { return (this->*_func_propose_)(prng, start); }
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

    for (const ushort &argID : _args_to_adapt_)
        (distr_ptr->set_arg)(argID, ((distr_ptr->get_arg)(argID)) * factor);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - show  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Base_Proposer::show_arg_names() const { return distr_ptr->show_arg_names(); }
const std::string Base_Proposer::show_arg_values() const { return distr_ptr->show_arg_values(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base Proposer - functional  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// alternative proposer methods
double Base_Proposer::_propose_undefined_(rnd_gen &prng, double start)
{
    tools::error("Attempt to call a method of an undefined proposer.");
    UNUSED(prng);
    UNUSED(start);
    return std::numeric_limits<double>::quiet_NaN();
}
double Base_Proposer::_propose_sum_to_previous_(rnd_gen &prng, double start) { return start + rand(prng); }
double Base_Proposer::_propose_times_previous_(rnd_gen &prng, double start) { return start * rand(prng); }
double Base_Proposer::_propose_draw_rand_(rnd_gen &prng, double start)
{
    UNUSED(start);
    return rand(prng);
}
// alternative logratios
double Base_Proposer::_logratio_undefined_(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.0;
}
double Base_Proposer::_logratio_pdf_(double xold, double xnew) const { return distr_ptr->logratio(xold, xnew); }
// FIXME : double check these
double Base_Proposer::_logratio_sum_to_zero_(double xold, double xnew) const { return _logratio_pdf_(xold - get_from_value(), xnew - get_from_value()); }
double Base_Proposer::_logratio_times_zero_(double xold, double xnew) const { return _logratio_pdf_(xold / get_from_value(), xnew / get_from_value()); }
double Base_Proposer::_logratio_symm_(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.0;
}
double Base_Proposer::_logratio_lognorm_symm_(double xold, double xnew) const { return log(xold) - log(xnew); }

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
        double rule_scale = 0.5,
        double rule_from = INF_DOUBLE,
        bool is_discr = false) : Base_Proposer(distribution("unif", dict_args(), dict_args(), is_discr),
                                               rule_from,
                                               {0, 1})
    {
        rule_scale = std::fabs(rule_scale);
        distr_ptr->assign_args({{"min", -rule_scale}, {"max", rule_scale}});

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
        double rule_scale = 1.0,
        double rule_from = INF_DOUBLE,
        bool is_discr = false) : Base_Proposer(distribution("norm", dict_args(), dict_args(), is_discr),
                                               rule_from,
                                               {1})
    {
        rule_scale = std::fabs(rule_scale);
        distr_ptr->assign_args({{"mean", 0.0}, {"sd", rule_scale}});

        _func_propose_ = &Normal_Proposer::_propose_sum_to_previous_;

        if (is_rule_start_from_zero())
            _logratio_method_ = &Normal_Proposer::_logratio_sum_to_zero_;
        else
            _logratio_method_ = &Normal_Proposer::_logratio_symm_;
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
        double rule_scale = 1.0,
        double rule_from = INF_DOUBLE,
        bool is_discr = false) : Base_Proposer(distribution("lnorm", dict_args(), dict_args(), is_discr),
                                               rule_from,
                                               {1})
    {
        rule_scale = std::fabs(rule_scale);
        distr_ptr->assign_args({{"meanlog", 0.0}, {"sdlog", rule_scale}});

        _func_propose_ = &Lognorm_Proposer::_propose_times_previous_;
        if (is_rule_start_from_zero())
            _logratio_method_ = &Lognorm_Proposer::_logratio_times_zero_;
        else
            _logratio_method_ = &Lognorm_Proposer::_logratio_lognorm_symm_;
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
std::shared_ptr<Base_Proposer> proposer(std::string rule_name,
                                        double rule_scale,
                                        double rule_from,
                                        bool is_discr)
{
    const std::vector<std::string> options = {"gibbs", "unif", "norm", "lnorm"};
    std::shared_ptr<Base_Proposer> output;

    if (tools::is_na(rule_name))
        output = std::make_shared<Base_Proposer>(std::make_shared<Base_Distribution>());

    if (rule_name == "unif")
        output = std::make_shared<Uniform_Proposer>(rule_scale, rule_from, is_discr);

    else if (rule_name == "norm")
        output = std::make_shared<Normal_Proposer>(rule_scale, rule_from, is_discr);

    else if (rule_name == "lnorm")
        output = std::make_shared<Lognorm_Proposer>(rule_scale, rule_from, is_discr);

    else
        throw std::invalid_argument(
            "Unknown proposer rule '" + rule_name + "'." + tools::find_closest(rule_name, options));

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Variable - getter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Variable::get_name() const { return name; }
double Base_Variable::get_value() const { return value; }
const double *Base_Variable::get_pointer() const { return &value; }
bool Base_Variable::get_infer() const { return infer; }
bool Base_Variable::is_from_zero() const { return (prop_ptr->is_rule_start_from_zero)(); }
double Base_Variable::get_from_value() const { return (prop_ptr->get_from_value)(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Variable - setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Base_Variable::Base_Variable(const std::string &name, double value, bool is_discr)
{
    this->name = name;
    this->value = value;
    this->is_discr = is_discr;
}

void Base_Variable::set_name(const std::string &x) { name = x; }
void Base_Variable::set_value(double x)
{
    buffer = value;
    value = x;
}
void Base_Variable::restore_value() { value = buffer; }

void Base_Variable::assign_distribution(const std::string &distr_name,
                                        dict_args distr_args,
                                        dict_args range)
{
    distr_ptr = distribution(distr_name, distr_args, range, is_discr);
}

void Base_Variable::assign_distribution(const std::string &distr_name,
                                        dict_args_ptr distr_args,
                                        dict_args range)
{
    distr_ptr = distribution(distr_name, distr_args, range, is_discr);
}

void Base_Variable::assign_proposer(const std::string &rule_name,
                                    double rule_scale,
                                    double rule_from)
{
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
            prop_ptr = proposer(rule_name, rule_scale, rule_from, is_discr);
            _func_propose_ = is_from_zero() ? &Base_Variable::_propose_from_zero_ : &Base_Variable::_propose_from_previous_;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Base_Variable - Metropolis-Hastings  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Variable::get_move_logratio() const { return (this->*_func_move_logratio_)(); }
bool Base_Variable::valid_move() const { return (this->*_func_valid_move_)(); }

//----- distr_ptr : prior/model methods -----//
void Base_Variable::rand(rnd_gen &prng) { set_value(distr_ptr->rand(prng)); }

double Base_Variable::distr_logpdf() const { return distr_ptr->logpdf(value); }
double Base_Variable::distr_logratio() const { return distr_ptr->logratio(buffer, value); }

//----- prop_ptr : propose methods -----//
void Base_Variable::propose(rnd_gen &prng) { return set_value((this->*_func_propose_)(prng)); }
double Base_Variable::propose_logratio() const { return (prop_ptr->logratio)(buffer, value); }

//----- prop_ptr : adaptive methods -----//
double Base_Variable::get_accept_rate() const { return prop_ptr->get_accept_rate(); }
void Base_Variable::update_accept_rate(bool is_accepted) { prop_ptr->update_accept_rate(is_accepted); }
void Base_Variable::reset_accept_rate() { prop_ptr->reset_accept_rate(); }
void Base_Variable::adapt_rule(double ratio) { (prop_ptr->adapt)(ratio); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base_Variable - show  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Base_Variable::show_distr_name() const { return distr_ptr->get_name(); }
const std::string Base_Variable::show_distr_arg_names() const { return distr_ptr->show_arg_names(); }
const std::string Base_Variable::show_distr_arg_values() const { return distr_ptr->show_arg_values(); }

const std::string Base_Variable::show_prop_name() const { return prop_ptr->get_name(); }
const std::string Base_Variable::show_prop_arg_names() const { return prop_ptr->show_arg_names(); }
const std::string Base_Variable::show_prop_arg_values() const { return prop_ptr->show_arg_values(); }

const std::string Base_Variable::show_range_min() const { return distr_ptr->show_range_min(); }
const std::string Base_Variable::show_range_max() const { return distr_ptr->show_range_max(); }

const std::string Base_Variable::show() const
{
    using namespace tools;
    std::ostringstream oss;
    oss << get_value() << C_DELIM << truth(infer);
    if (infer)
        oss << C_DELIM << show_distr_name() << C_DELIM << show_distr_arg_names() << C_DELIM << show_distr_arg_values()
            << C_DELIM << show_prop_name() << C_DELIM << show_prop_arg_names() << C_DELIM << show_prop_arg_values()
            << C_DELIM << (tools::is_na(show_range_min()) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(show_range_min()))
            << C_DELIM << (tools::is_na(show_range_max()) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(show_range_max()));

    return oss.str();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base_Variable - functional  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

bool Base_Variable::_valid_move_default_() const { return distr_ptr->valid(value); }
bool Base_Variable::_valid_move_gibbs_() const { return true; }

double Base_Variable::_move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
double Base_Variable::_move_logratio_gibbs_() const { return 0.0; }

double Base_Variable::_propose_error_(rnd_gen &prng)
{
    UNUSED(prng);
    tools::error("Call to proposer of immutable variable '" + get_name() + "'.");
    return 0.;
}
double Base_Variable::_propose_gibbs_sample_(rnd_gen &prng) { return (prop_ptr->propose)(prng, 0); }
double Base_Variable::_propose_from_previous_(rnd_gen &prng) { return (prop_ptr->propose)(prng, get_value()); }
double Base_Variable::_propose_from_zero_(rnd_gen &prng) { return (prop_ptr->propose)(prng, get_from_value()); }
