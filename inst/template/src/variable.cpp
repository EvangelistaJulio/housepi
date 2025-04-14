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
                   bool is_discrete) : _distr_ptr_(distr)
{
    // set to discrete
    _distr_ptr_->init_reshaping(is_discrete);

    // assign user rule scale
    if (!tools::is_na(rule_scale))
        _user_rule_scale_ = std::fabs(rule_scale);

    // by default proposal is from previous state
    _is_start_from_zero_ = false;

    // proposal alwats starting by given value
    if (!tools::is_na(rule_start_value))
    {
        _is_start_from_zero_ = true;
        set_rule_start_value(rule_start_value);
    }

    // load which are the parameter to be adapted
    n_args_to_adapt = args_to_adapt.size();
    _args_to_adapt_ = new ushort[n_args_to_adapt];

    for (ushort c = 0; c < n_args_to_adapt; ++c)
        _args_to_adapt_[c] = args_to_adapt.at(c);
};

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Proposer - getter  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Proposer::show_name() const { return this->name; }
double Proposer::get_arg(ushort argID) const { return (_distr_ptr_->get_arg)(argID); }
double Proposer::get_from_value() const { return from_value; }
bool Proposer::is_rule_start_from_zero() const { return _is_start_from_zero_; }
double Proposer::propose(rnd_gen &prng, double start) { return (this->*_func_propose_)(prng, start); }
double Proposer::logratio(double xold, double xnew) const { return (this->*_func_logratio_)(xold, xnew); }
double Proposer::rand(rnd_gen &prng) { return (_distr_ptr_->rand)(prng); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Proposer - setter  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Proposer::user_reset() {};
void Proposer::set_rule_start_value(double x) { from_value = x; }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base Proposer - adaptation  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Proposer::get_accept_rate() const { return accept_rate; }
void Proposer::update_accept_rate(bool is_accepted)
{
    accept_rate = (double)(accept_rate * n_proposed_moves + is_accepted) / (n_proposed_moves + 1);
    ++n_proposed_moves;
}
void Proposer::reset_accept_rate()
{
    accept_rate = 0.;
    n_proposed_moves = 0;
};

/* Increase spread if acceptation rate is high and viceversa */
void Proposer::adapt(double ratio)
{
    double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);
    // loop over all arguments to adapt
    for (ushort c = 0; c < n_args_to_adapt; ++c)
    {
        ushort arg_id = _args_to_adapt_[c];
        (_distr_ptr_->set_arg)(arg_id,
                               (_distr_ptr_->get_arg)(arg_id)*factor);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Show - adaptation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Proposer::show_arg_names() const { return _distr_ptr_->show_arg_names(); }
const std::string Proposer::show_arg_values() const { return _distr_ptr_->show_arg_values(); }
const std::string Proposer::show_scale() const
{
    if (n_args_to_adapt == 0)
        return NULL_STRING;
    else
        return std::to_string(std::fabs((_distr_ptr_->get_arg)(_args_to_adapt_[0])));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base Proposer - functional  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// alternative proposer methods
double Proposer::_propose_undefined_(rnd_gen &prng, double start)
{
    UNUSED(prng);
    UNUSED(start);
    throw std::invalid_argument(
        "Attempt to call a method of an undefined proposer.");
    return std::numeric_limits<double>::quiet_NaN();
}
double Proposer::_propose_sum_to_previous_(rnd_gen &prng, double start) { return start + rand(prng); }
double Proposer::_propose_times_previous_(rnd_gen &prng, double start) { return start * rand(prng); }
double Proposer::_propose_draw_rand_(rnd_gen &prng, double start)
{
    UNUSED(start);
    return rand(prng);
}
// alternative logratios
double Proposer::_logratio_undefined_(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.;
}
double Proposer::_logratio_pdf_(double xold, double xnew) const { return _distr_ptr_->logratio(xold, xnew); }
double Proposer::_logratio_sum_to_zero_(double xold, double xnew) const { return _logratio_pdf_(xold - from_value, xnew - from_value); }
double Proposer::_logratio_times_zero_(double xold, double xnew) const { return _logratio_pdf_(xold / from_value, xnew / from_value); }
double Proposer::_logratio_symm_(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.;
}
double Proposer::_logratio_lognorm_symm_(double xold, double xnew) const { return log(xold) - log(xnew); }

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
        _distr_ptr_ = distr;

        name = "draw";
        // proposer functional method
        _func_propose_ = &Draw_Proposer::_propose_draw_rand_;
        // logratio functional method
        _func_logratio_ = &Draw_Proposer::_logratio_pdf_;
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
        bool is_discrete = false) : Proposer(distribution_switchboard("unif"),
                                             rule_scale,
                                             rule_from,
                                             {0, 1},
                                             is_discrete)
    {
        name = _distr_ptr_->show_name();
        // initialize
        user_reset();
        // proposer functional method
        _func_propose_ = &Uniform_Proposer::_propose_sum_to_previous_;
        // logratio functional method
        _func_logratio_ = &Uniform_Proposer::_logratio_symm_;
    }
    ~Uniform_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr_->assign_scalar_arg("min", -_user_rule_scale_);
        _distr_ptr_->assign_scalar_arg("max", _user_rule_scale_);
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
        bool is_discrete = false) : Proposer(distribution_switchboard("norm"),
                                             rule_scale,
                                             rule_from,
                                             {1},
                                             is_discrete)
    {
        name = _distr_ptr_->show_name();
        // initialize
        user_reset();
        // proposer functional method
        _func_propose_ = &Normal_Proposer::_propose_sum_to_previous_;
        // logratio functional method
        if (is_rule_start_from_zero())
            _func_logratio_ = &Normal_Proposer::_logratio_sum_to_zero_;
        else
            _func_logratio_ = &Normal_Proposer::_logratio_symm_;
    }
    ~Normal_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr_->assign_scalar_arg("mean", 0.);
        _distr_ptr_->assign_scalar_arg("sd", _user_rule_scale_);
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
        bool is_discrete = false) : Proposer(distribution_switchboard("lnorm"),
                                             rule_scale,
                                             rule_from,
                                             {1},
                                             is_discrete)
    {
        name = _distr_ptr_->show_name();
        // initialize
        user_reset();
        // proposer functional method
        _func_propose_ = &Lognorm_Proposer::_propose_times_previous_;
        // logratio functional method
        if (is_rule_start_from_zero())
            _func_logratio_ = &Lognorm_Proposer::_logratio_times_zero_;
        else
            _func_logratio_ = &Lognorm_Proposer::_logratio_lognorm_symm_;
    }
    ~Lognorm_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr_->assign_scalar_arg("meanlog", 0.);
        _distr_ptr_->assign_scalar_arg("sdlog", _user_rule_scale_);
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Switchboard to select the custom metropolis proposer rule. */
std::shared_ptr<Proposer> proposer_switchboard(const std::string &rule_name = NULL_STRING,
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
bool Variable::_valid_move_default_() const { return _distr_ptr_->valid(value); }
bool Variable::_valid_move_draw_() const { return true; }

//===== Move logratio : functional methods =====//
double Variable::_move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
double Variable::_move_logratio_draw_() const { return 0.; }

//===== Propose : functional methods =====//
double Variable::_propose_error_(rnd_gen &prng)
{
    UNUSED(prng);
    throw std::invalid_argument(
        "Call to proposer of constant variable \"" + name + "\".");
    return 0.;
}
double Variable::_propose_draw_sample_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, 0); }
double Variable::_propose_from_previous_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, value); }
double Variable::_propose_from_zero_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, _prop_ptr_->get_from_value()); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Variable - METROPOLIS-HASTINGS  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store the buffer and change value. */
void Variable::set_value(double x)
{
    buffer = value;
    value = x;
}

/* Restore buffered value. */
void Variable::restore_value() { value = buffer; }

double Variable::get_move_logratio() const { return (this->*_func_move_logratio_)(); }
bool Variable::valid_move() const { return (this->*_func_valid_move_)(); }

//----- distr_ptr : prior/model methods -----//
void Variable::rand(rnd_gen &prng) { set_value(_distr_ptr_->rand(prng)); }

double Variable::distr_logpdf() const { return _distr_ptr_->logpdf(value); }
double Variable::distr_logratio() const { return _distr_ptr_->logratio(buffer, value); }

//----- prop_ptr : propose methods -----//
void Variable::propose(rnd_gen &prng) { set_value((this->*_func_propose_)(prng)); }
double Variable::propose_logratio() const { return (_prop_ptr_->logratio)(buffer, value); }

//----- prop_ptr : adaptive methods -----//
double Variable::get_accept_rate() const { return _prop_ptr_->get_accept_rate(); }
void Variable::update_accept_rate(bool is_accepted) { _prop_ptr_->update_accept_rate(is_accepted); }
void Variable::reset_accept_rate() { _prop_ptr_->reset_accept_rate(); }
void Variable::adapt_rule(double ratio) { (_prop_ptr_->adapt)(ratio); }

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
                   bool is_discrete) : _distr_ptr_(std::make_shared<Distribution>()),
                                       _prop_ptr_(std::make_shared<Proposer>())
{
    init(name, value, is_discrete);
}
Variable::~Variable() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable - GETTER  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Variable::get() const { return value; }
const double *Variable::get_pointer() const { return &value; }
bool Variable::is_constant() const { return _is_constant_; }
bool Variable::is_discrete() const { return _distr_ptr_->is_discrete(); }
bool Variable::is_from_zero() const { return (_prop_ptr_->is_rule_start_from_zero)(); }
double Variable::get_from_value() const { return (_prop_ptr_->get_from_value)(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable - SETTER  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Restore user defined variable initialization. */
void Variable::init(const std::string &name, double value, bool is_discrete)
{
    this->name = name;
    this->_user_value_ = value;
    this->value = value;
    this->_is_discrete_ = is_discrete;
}

/* Restore user defined variable initialization. */
void Variable::user_reset()
{
    this->value = _user_value_;
    this->_prop_ptr_->user_reset();
}

/* Assign distribution with constant arguments. */
void Variable::assign_distribution(const std::string &distr_name)
{
    _distr_ptr_ = distribution_switchboard(distr_name);
}

void Variable::assign_distr_scalar_arg(const std::string &arg_name, double x)
{
    _distr_ptr_->assign_scalar_arg(arg_name, x);
}
void Variable::assign_distr_pointer_arg(const std::string &arg_name, const double *y)
{
    _distr_ptr_->assign_pointer_arg(arg_name, y);
}

void Variable::assign_distr_scalar_range_min(double x)
{
    _distr_ptr_->assign_scalar_range_min(x);
}
void Variable::assign_distr_pointer_range_min(const double *y)
{
    _distr_ptr_->assign_pointer_range_min(y);
}

void Variable::assign_distr_scalar_range_max(double x)
{
    _distr_ptr_->assign_scalar_range_max(x);
}
void Variable::assign_distr_pointer_range_max(const double *y)
{
    _distr_ptr_->assign_pointer_range_max(y);
}

void Variable::init_distr_reshaping()
{
    _distr_ptr_->init_reshaping(_is_discrete_);
}

/* Assign proposer. */
void Variable::assign_proposer(const std::string &rule_name,
                               double rule_scale,
                               double rule_from)
{
    if (tools::is_na(rule_name))
    {
        // if no rule is defined constant variable
        _is_constant_ = true;
        _func_propose_ = &Variable::_propose_error_;
        _prop_ptr_ = proposer_switchboard();
    }
    else
    {
        // variable to sample
        _is_constant_ = false;

        // define proposer rule with functional programming
        if (rule_name == "draw")
        {
            //----- Draw proposer rule -----//
            // assign the pointer to the proposer
            _prop_ptr_ = std::make_shared<Draw_Proposer>(_distr_ptr_);

            // sampler
            _func_propose_ = &Variable::_propose_draw_sample_;
            // logratio
            _func_move_logratio_ = &Variable::_move_logratio_draw_;
            // validation
            _func_valid_move_ = &Variable::_valid_move_draw_;
        }
        else
        {
            //----- Other proposer rules -----//
            // value from which to start proposal
            if ((!tools::is_na(rule_from)) && (!_distr_ptr_->valid(rule_from)))
                tools::warning(
                    "Proposer rule from value " + std::to_string(rule_from) + " out of range" + " {min:" + show_min() + ", max:" + show_max() + "}.");

            // assign the pointer to the proposer
            _prop_ptr_ = proposer_switchboard(rule_name, rule_scale, rule_from, is_discrete());

            // raise a warning if "lnorm" rule when varaible can be negative
            if ((_distr_ptr_->get_min() < 0) && (rule_name == "lnorm"))
                tools::warning(
                    "Lognormal proposer rule for non-positive variable \"" + name + "\" .");

            // sampler
            if (is_from_zero())
                _func_propose_ = &Variable::_propose_from_zero_;
            else
                _func_propose_ = &Variable::_propose_from_previous_;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable - SHOW  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Variable::show_name() const { return name; }
const std::string Variable::show_distr_name() const { return _distr_ptr_->show_name(); }
const std::string Variable::show_distr_arg_names() const { return _distr_ptr_->show_arg_names(); }
const std::string Variable::show_distr_arg_values() const { return _distr_ptr_->show_arg_values(); }

const std::string Variable::show_prop_name() const { return _prop_ptr_->show_name(); }
const std::string Variable::show_prop_arg_names() const { return _prop_ptr_->show_arg_names(); }
const std::string Variable::show_prop_arg_values() const { return _prop_ptr_->show_arg_values(); }
const std::string Variable::show_prop_scale() const { return _prop_ptr_->show_scale(); }

const std::string Variable::show_min() const { return _distr_ptr_->show_min(); }
const std::string Variable::show_max() const { return _distr_ptr_->show_max(); }
