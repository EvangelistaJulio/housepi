/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * variable.cpp
 *
 */

#include "variable.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Proposer : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Proposer::Proposer(std::shared_ptr<Distribution> distr,
                   double rule_scale,
                   double rule_start_value,
                   const std::vector<ushort> &args_to_adapt,
                   bool is_discrete) : _distr_ptr(distr)
{
    // default proposer name derives from distribution
    _proposer_name = _distr_ptr->show_name();

    // set to discrete
    _distr_ptr->init_reshaping(is_discrete);

    // assign user rule scale
    if (!tools::is_na(rule_scale))
        _user_rule_scale = std::fabs(rule_scale);

    // by default proposal is from previous state
    _is_start_from_zero = false;

    // proposal alwats starting by given value
    if (!tools::is_na(rule_start_value))
    {
        _is_start_from_zero = true;
        set_rule_start_value(rule_start_value);
    }

    // load which are the parameter to be adapted
    _n_args_to_adapt = args_to_adapt.size();
    _args_to_adapt_list = new ushort[_n_args_to_adapt];

    for (ushort c = 0; c < _n_args_to_adapt; ++c)
        _args_to_adapt_list[c] = args_to_adapt.at(c);
};

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Proposer - getter  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Proposer::show_name() const { return this->_proposer_name; }
double Proposer::get_arg(ushort argID) const { return (_distr_ptr->get_arg)(argID); }
double Proposer::get_from_value() const { return _from_value; }
bool Proposer::is_rule_start_from_zero() const { return _is_start_from_zero; }
double Proposer::propose(rnd_gen &prng, double start) { return (this->*_func_propose)(prng, start); }
double Proposer::logratio(double xold, double xnew) const { return (this->*_func_logratio)(xold, xnew); }
double Proposer::rand(rnd_gen &prng) { return (_distr_ptr->rand)(prng); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Proposer - setter  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Proposer::user_reset() {};
void Proposer::set_rule_start_value(double x) { _from_value = x; }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base Proposer - adaptation  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Proposer::get_accept_rate() const { return _accept_rate; }
void Proposer::update_accept_rate(bool is_accepted)
{
    _accept_rate = (double)(_accept_rate * _n_proposed_moves + is_accepted) / (_n_proposed_moves + 1);
    ++_n_proposed_moves;
}
void Proposer::reset_accept_rate()
{
    _accept_rate = 0.;
    _n_proposed_moves = 0;
};

/* Increase spread if acceptation rate is high and viceversa */
void Proposer::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);
    // loop over all arguments to adapt
    for (ushort c = 0; c < _n_args_to_adapt; ++c)
    {
        ushort arg_id = _args_to_adapt_list[c];
        (_distr_ptr->set_arg)(arg_id, (_distr_ptr->get_arg)(arg_id)*factor);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Show - adaptation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Proposer::show_arg_names() const { return _distr_ptr->show_arg_names(); }
const std::string Proposer::show_arg_values() const { return _distr_ptr->show_arg_values(); }
const std::string Proposer::show_scale() const
{
    if (_n_args_to_adapt == 0)
        return NULL_STRING;
    else
        return std::to_string(std::fabs((_distr_ptr->get_arg)(_args_to_adapt_list[0])));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base Proposer - functional  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// alternative proposer methods
double Proposer::__propose_undefined(rnd_gen &prng, double start)
{
    UNUSED(prng);
    UNUSED(start);
    throw std::invalid_argument(
        "Attempt to call a method of an undefined proposer.");
    return std::numeric_limits<double>::quiet_NaN();
}
double Proposer::__propose_sum_to_previous(rnd_gen &prng, double start) { return start + rand(prng); }
double Proposer::__propose_times_previous(rnd_gen &prng, double start) { return start * rand(prng); }
double Proposer::__propose_draw_rand(rnd_gen &prng, double start)
{
    UNUSED(start);
    return rand(prng);
}
// alternative logratios
double Proposer::__logratio_undefined(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.;
}
double Proposer::__logratio_pdf(double xold, double xnew) const { return _distr_ptr->logratio(xold, xnew); }
double Proposer::__logratio_sum_to_zero(double xold, double xnew) const { return __logratio_pdf(xold - _from_value, xnew - _from_value); }
double Proposer::__logratio_times_zero(double xold, double xnew) const { return __logratio_pdf(xold / _from_value, xnew / _from_value); }
double Proposer::__logratio_symm(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.;
}
double Proposer::__logratio_lognorm_symm(double xold, double xnew) const { return std::log(xold / xnew); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Draw_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Draw_Proposer : public Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Draw_Proposer(std::shared_ptr<Distribution> distr)
    {
        _distr_ptr = distr;

        _proposer_name = "draw";
        // proposer functional method
        _func_propose = &Draw_Proposer::__propose_draw_rand;
        // logratio functional method
        _func_logratio = &Draw_Proposer::__logratio_pdf;
    }
    ~Draw_Proposer() = default;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Uniform_Proposer(
        double rule_scale = NULL_DOUBLE,
        double rule_from = NULL_DOUBLE,
        bool is_discrete = false) : Proposer(distribution_switch("unif"),
                                             rule_scale,
                                             rule_from,
                                             {0, 1},
                                             is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        _func_propose = &Uniform_Proposer::__propose_sum_to_previous;
        // logratio functional method
        _func_logratio = &Uniform_Proposer::__logratio_symm;
    }
    ~Uniform_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr->assign_scalar_arg("min", -_user_rule_scale);
        _distr_ptr->assign_scalar_arg("max", _user_rule_scale);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Normal_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Normal_Proposer(
        double rule_scale = NULL_DOUBLE,
        double rule_from = NULL_DOUBLE,
        bool is_discrete = false) : Proposer(distribution_switch("norm"),
                                             rule_scale,
                                             rule_from,
                                             {1},
                                             is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        _func_propose = &Normal_Proposer::__propose_sum_to_previous;
        // logratio functional method
        if (is_rule_start_from_zero())
            _func_logratio = &Normal_Proposer::__logratio_sum_to_zero;
        else
            _func_logratio = &Normal_Proposer::__logratio_symm;
    }
    ~Normal_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr->assign_scalar_arg("mean", 0.);
        _distr_ptr->assign_scalar_arg("sd", _user_rule_scale);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Lognorm_Proposer(
        double rule_scale = 1.,
        double rule_from = NULL_DOUBLE,
        bool is_discrete = false) : Proposer(distribution_switch("lnorm"),
                                             rule_scale,
                                             rule_from,
                                             {1},
                                             is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        _func_propose = &Lognorm_Proposer::__propose_times_previous;
        // logratio functional method
        if (is_rule_start_from_zero())
            _func_logratio = &Lognorm_Proposer::__logratio_times_zero;
        else
            _func_logratio = &Lognorm_Proposer::__logratio_lognorm_symm;
    }
    ~Lognorm_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr->assign_scalar_arg("meanlog", 0.);
        _distr_ptr->assign_scalar_arg("sdlog", _user_rule_scale);
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Switchboard to select the custom metropolis proposer rule. */
std::shared_ptr<Proposer> proposer_switch(const std::string &rule_name = NULL_STRING,
                                          double rule_scale = NULL_DOUBLE,
                                          double rule_from = NULL_DOUBLE,
                                          bool is_discrete = false)
{
    const str_vec options = {"draw", "unif", "norm", "lnorm"};

    if (tools::is_na(rule_name))
        return std::make_shared<Proposer>(std::make_shared<Distribution>());

    else if (rule_name == "unif")
        return std::make_shared<Uniform_Proposer>(rule_scale, rule_from, is_discrete);

    else if (rule_name == "norm")
        return std::make_shared<Normal_Proposer>(rule_scale, rule_from, is_discrete);

    else if (rule_name == "lnorm")
        return std::make_shared<Lognorm_Proposer>(rule_scale, rule_from, is_discrete);

    else
        throw std::invalid_argument(
            "Unknown proposer rule \"" + rule_name + "\"." + tools::find_closest(rule_name, options));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Variable - functional  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== Valid move : functional methods =====//
bool Variable::__valid_move_default() const { return _distr_ptr->valid(_value); }
bool Variable::__valid_move_draw() const { return true; }

//===== Move logratio : functional methods =====//
double Variable::__move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
double Variable::__move_logratio_draw() const { return 0.; }

//===== Propose : functional methods =====//
double Variable::__propose_error(rnd_gen &prng)
{
    UNUSED(prng);
    throw std::invalid_argument(
        "Call to proposer of constant variable \"" + _variable_name + "\".");
    return 0.;
}
double Variable::__propose_draw_sample(rnd_gen &prng) { return (_prop_ptr->propose)(prng, 0); }
double Variable::__propose_from_previous(rnd_gen &prng) { return (_prop_ptr->propose)(prng, _value); }
double Variable::__propose_from_zero(rnd_gen &prng) { return (_prop_ptr->propose)(prng, _prop_ptr->get_from_value()); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Variable - METROPOLIS-HASTINGS  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store the buffer and change value. */
void Variable::set_value(double x)
{
    _buffer_value = _value;
    _value = x;
}

/* Restore buffered value. */
void Variable::restore_value() { _value = _buffer_value; }

double Variable::get_move_logratio() const { return (this->*_func_move_logratio)(); }
bool Variable::valid_move() const { return (this->*_func_valid_move)(); }

//----- distr_ptr : prior/model methods -----//
void Variable::rand(rnd_gen &prng) { set_value(_distr_ptr->rand(prng)); }

double Variable::distr_logpdf() const { return _distr_ptr->logpdf(_value); }
double Variable::distr_logratio() const { return _distr_ptr->logratio(_buffer_value, _value); }

//----- prop_ptr : propose methods -----//
void Variable::propose(rnd_gen &prng) { set_value((this->*_func_propose)(prng)); }
double Variable::propose_logratio() const { return (_prop_ptr->logratio)(_buffer_value, _value); }

//----- prop_ptr : adaptive methods -----//
double Variable::get_accept_rate() const { return _prop_ptr->get_accept_rate(); }
void Variable::update_accept_rate(bool is_accepted) { _prop_ptr->update_accept_rate(is_accepted); }
void Variable::reset_accept_rate() { _prop_ptr->reset_accept_rate(); }
void Variable::adapt_rule(double ratio) { (_prop_ptr->adapt)(ratio); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Variable - CONSTRUCTOR  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Variable::Variable(const std::string &name,
                   double value,
                   bool is_discrete) : _distr_ptr(std::make_shared<Distribution>()),
                                       _prop_ptr(std::make_shared<Proposer>())
{
    init(name, value, is_discrete);
}
Variable::~Variable() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable - GETTER  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Variable::get() const { return _value; }
const double *Variable::get_pointer() const { return &_value; }
bool Variable::is_constant() const { return _is_constant; }
bool Variable::is_discrete() const { return _distr_ptr->is_discrete(); }
bool Variable::is_from_zero() const { return (_prop_ptr->is_rule_start_from_zero)(); }
double Variable::get_from_value() const { return (_prop_ptr->get_from_value)(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable - SETTER  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Restore user defined variable initialization. */
void Variable::init(const std::string &name, double value, bool is_discrete)
{
    _variable_name = name;
    _user_value = value;
    _value = _user_value;
    _is_discrete = is_discrete;
}

/* Restore user defined variable initialization. */
void Variable::user_reset()
{
    _value = _user_value;
    _prop_ptr->user_reset();
}

/* Assign distribution with constant arguments. */
void Variable::assign_distribution(const std::string &distr_name)
{
    _distr_ptr = distribution_switch(distr_name);
}

void Variable::assign_distr_scalar_arg(const std::string &arg_name, double x)
{
    _distr_ptr->assign_scalar_arg(arg_name, x);
}
void Variable::assign_distr_pointer_arg(const std::string &arg_name, const double *y)
{
    _distr_ptr->assign_pointer_arg(arg_name, y);
}

void Variable::assign_distr_scalar_range_min(double x)
{
    _distr_ptr->assign_scalar_range_min(x);
}
void Variable::assign_distr_pointer_range_min(const double *y)
{
    _distr_ptr->assign_pointer_range_min(y);
}

void Variable::assign_distr_scalar_range_max(double x)
{
    _distr_ptr->assign_scalar_range_max(x);
}
void Variable::assign_distr_pointer_range_max(const double *y)
{
    _distr_ptr->assign_pointer_range_max(y);
}

void Variable::init_distr_reshaping()
{
    _distr_ptr->init_reshaping(_is_discrete);
}

/* Assign proposer. */
void Variable::assign_proposer(const std::string &rule_name,
                               double rule_scale,
                               double rule_from)
{
    if (tools::is_na(rule_name))
    {
        // if no rule is defined constant variable
        _is_constant = true;
        _func_propose = &Variable::__propose_error;
        _prop_ptr = proposer_switch();
    }
    else
    {
        // variable to sample
        _is_constant = false;

        // define proposer rule with functional programming
        if (rule_name == "draw")
        {
            //----- Draw proposer rule -----//
            // assign the pointer to the proposer
            _prop_ptr = std::make_shared<Draw_Proposer>(_distr_ptr);

            // sampler
            _func_propose = &Variable::__propose_draw_sample;
            // logratio
            _func_move_logratio = &Variable::__move_logratio_draw;
            // validation
            _func_valid_move = &Variable::__valid_move_draw;
        }
        else
        {
            //----- Other proposer rules -----//
            // value from which to start proposal
            if ((!tools::is_na(rule_from)) && (!_distr_ptr->valid(rule_from)))
                tools::warning(
                    "Proposer rule from value " + std::to_string(rule_from) + " out of range" + " {min:" + show_min() + ", max:" + show_max() + "}.");

            // assign the pointer to the proposer
            _prop_ptr = proposer_switch(rule_name, rule_scale, rule_from, is_discrete());

            // raise a warning if "lnorm" rule when varaible can be negative
            if ((_distr_ptr->get_min() < 0) && (rule_name == "lnorm"))
                tools::warning(
                    "Lognormal proposer rule for non-positive variable \"" + _variable_name + "\" .");

            // sampler
            if (is_from_zero())
                _func_propose = &Variable::__propose_from_zero;
            else
                _func_propose = &Variable::__propose_from_previous;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable - SHOW  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Variable::show_name() const { return _variable_name; }
const std::string Variable::show_distr_name() const { return _distr_ptr->show_name(); }
const std::string Variable::show_distr_arg_names() const { return _distr_ptr->show_arg_names(); }
const std::string Variable::show_distr_arg_values() const { return _distr_ptr->show_arg_values(); }

const std::string Variable::show_prop_name() const { return _prop_ptr->show_name(); }
const std::string Variable::show_prop_arg_names() const { return _prop_ptr->show_arg_names(); }
const std::string Variable::show_prop_arg_values() const { return _prop_ptr->show_arg_values(); }
const std::string Variable::show_prop_scale() const { return _prop_ptr->show_scale(); }

const std::string Variable::show_min() const { return _distr_ptr->show_min(); }
const std::string Variable::show_max() const { return _distr_ptr->show_max(); }
