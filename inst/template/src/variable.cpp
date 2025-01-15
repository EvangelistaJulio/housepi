/*
 * variable.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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
                             double rule_scale,
                             double rule_start_value,
                             const std::vector<ushort> &args_to_adapt,
                             bool is_discrete) : _distr_ptr_(distr)
{
    // set to discrete
    _distr_ptr_->assign_reshaping(is_discrete,
                                  NULL_DOUBLE,
                                  NULL_DOUBLE);

    // assign user rule scale
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
    _args_to_adapt_ = args_to_adapt;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - getter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Proposer::get_name() const { return this->name; }
double Base_Proposer::get_arg(ushort argID) const { return (_distr_ptr_->get_arg)(argID); }
double Base_Proposer::get_from_value() const { return from_value; }
bool Base_Proposer::is_rule_start_from_zero() const { return _is_start_from_zero_; }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Base_Proposer::user_reset() {}
void Base_Proposer::set_rule_start_value(double x) { from_value = x; }
double Base_Proposer::rand(rnd_gen &prng) { return (_distr_ptr_->rand)(prng); }
double Base_Proposer::propose(rnd_gen &prng, double start) { return (this->*_func_propose_)(prng, start); }
double Base_Proposer::logratio(double xold, double xnew) const { return (this->*_func_logratio_)(xold, xnew); }

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Proposer - adaptation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Proposer::get_accept_rate() const { return accept_rate; }

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
        (_distr_ptr_->set_arg)(argID, (_distr_ptr_->get_arg)(argID)*factor);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Proposer - show  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Base_Proposer::show_arg_names() const { return _distr_ptr_->show_arg_names(); }
const std::string Base_Proposer::show_arg_values() const { return _distr_ptr_->show_arg_values(); }
const std::string Base_Proposer::show_scale() const
{
    if (_args_to_adapt_.empty())
        return NULL_STRING;
    else
        return std::to_string(std::fabs((_distr_ptr_->get_arg)(_args_to_adapt_[0])));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base Proposer - functional  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

// alternative proposer methods
double Base_Proposer::_propose_undefined_(rnd_gen &prng, double start)
{
    UNUSED(prng);
    UNUSED(start);
    throw std::invalid_argument(
        "Attempt to call a method of an undefined proposer.");
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
double Base_Proposer::_logratio_pdf_(double xold, double xnew) const { return _distr_ptr_->logratio(xold, xnew); }
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
///////////////////////////////////  Gibbs_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gibbs_Proposer : public Base_Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Gibbs_Proposer(std::shared_ptr<Base_Distribution> distr) : Base_Proposer(distr)
    {
        name = "gibbs";
        // method for proposer
        _func_propose_ = &Gibbs_Proposer::_propose_draw_rand_;
        // method for logratio
        _func_logratio_ = &Gibbs_Proposer::_logratio_pdf_;
    }
    ~Gibbs_Proposer() = default;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Base_Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Uniform_Proposer(
        double rule_scale = 0.5,
        double rule_from = NULL_DOUBLE,
        bool is_discrete = false) : Base_Proposer(distribution_switchboard("unif"),
                                                  rule_scale,
                                                  rule_from,
                                                  {0, 1},
                                                  is_discrete)
    {
        name = _distr_ptr_->get_name();
        // initialize
        user_reset();
        // method for proposer
        _func_propose_ = &Uniform_Proposer::_propose_sum_to_previous_;
        // method for logratio
        _func_logratio_ = &Uniform_Proposer::_logratio_symm_;
    }
    ~Uniform_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override
    {
        _distr_ptr_->assign_scalar_arg("min", -_user_rule_scale_);
        _distr_ptr_->assign_scalar_arg("max", _user_rule_scale_);
        _distr_ptr_->set_range(); 
        // FIXME : if user_rule_scale changes, the range doesn't!
        // FIXME: time for range to be function arguments...
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Normal_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Base_Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Normal_Proposer(
        double rule_scale = 1.0,
        double rule_from = NULL_DOUBLE,
        bool is_discrete = false) : Base_Proposer(distribution_switchboard("norm"),
                                                  rule_scale,
                                                  rule_from,
                                                  {1},
                                                  is_discrete)
    {
        name = _distr_ptr_->get_name();
        // initialize
        user_reset();
        // method for proposer
        _func_propose_ = &Normal_Proposer::_propose_sum_to_previous_;
        // method for logratio
        if (is_rule_start_from_zero())
            _func_logratio_ = &Normal_Proposer::_logratio_sum_to_zero_;
        else
            _func_logratio_ = &Normal_Proposer::_logratio_symm_;
    }
    ~Normal_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override
    {
        _distr_ptr_->assign_scalar_arg("mean", 0.0);
        _distr_ptr_->assign_scalar_arg("sd", _user_rule_scale_);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Base_Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Lognorm_Proposer(
        double rule_scale = 1.0,
        double rule_from = NULL_DOUBLE,
        bool is_discrete = false) : Base_Proposer(distribution_switchboard("lnorm"),
                                                  rule_scale,
                                                  rule_from,
                                                  {1},
                                                  is_discrete)
    {
        name = _distr_ptr_->get_name();
        // initialize
        user_reset();
        // method for proposer
        _func_propose_ = &Lognorm_Proposer::_propose_times_previous_;
        // method for logratio
        if (is_rule_start_from_zero())
            _func_logratio_ = &Lognorm_Proposer::_logratio_times_zero_;
        else
            _func_logratio_ = &Lognorm_Proposer::_logratio_lognorm_symm_;
    }
    ~Lognorm_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override
    {
        _distr_ptr_->assign_scalar_arg("meanlog", 0.0);
        _distr_ptr_->assign_scalar_arg("sdlog", _user_rule_scale_);
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

const str_vec AVAILABLE_PROPOSER_RULE = {
    "gibbs", "unif", "norm", "lnorm"};

/* Switchboard to select the custom metropolis proposer rule. */
std::shared_ptr<Base_Proposer> proposer_switchboard(const std::string &rule_name = NULL_STRING,
                                                    double rule_scale = NULL_DOUBLE,
                                                    double rule_from = NULL_DOUBLE,
                                                    bool is_discrete = false)
{
    const str_vec &options = AVAILABLE_PROPOSER_RULE;

    if (!tools::is_na(rule_name))
    {
        if (rule_name == "unif")
            return std::make_shared<Uniform_Proposer>(rule_scale, rule_from, is_discrete);

        else if (rule_name == "norm")
            return std::make_shared<Normal_Proposer>(rule_scale, rule_from, is_discrete);

        else if (rule_name == "lnorm")
            return std::make_shared<Lognorm_Proposer>(rule_scale, rule_from, is_discrete);

        else
            throw std::invalid_argument(
                "Unknown proposer rule '" + rule_name + "'." + tools::find_closest(rule_name, options));
    }

    // otherwise default
    return std::make_shared<Base_Proposer>(std::make_shared<Base_Distribution>());
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base_Variable - CONSTRUCTOR  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Base_Variable::Base_Variable(const std::string &name, double value, bool is_discrete)
{
    init(name, value, is_discrete);
    assign_distribution();
    assign_proposer();
}

Base_Variable::~Base_Variable() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Variable - GETTER  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Variable::get_name() const { return name; }
double Base_Variable::get_value() const { return value; }
const double *Base_Variable::get_pointer() const { return &value; }
bool Base_Variable::is_constant() const { return _is_constant_; }
bool Base_Variable::is_discrete() const { return _is_discrete_; }
bool Base_Variable::is_from_zero() const { return (_prop_ptr_->is_rule_start_from_zero)(); }
double Base_Variable::get_from_value() const { return (_prop_ptr_->get_from_value)(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Variable - SETTER  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Restore user defined variable initialization. */
void Base_Variable::init(const std::string &name, double value, bool is_discrete)
{
    this->name = name;
    this->_user_value_ = value;
    this->value = value;
    this->_is_discrete_ = is_discrete;
}

/* Restore user defined variable initialization. */
void Base_Variable::user_reset()
{
    value = _user_value_;
    _prop_ptr_->user_reset();
}

/* Assign distribution with constant arguments. */
void Base_Variable::assign_distribution(const std::string &distr_name)
{
    _distr_ptr_ = distribution_switchboard(distr_name);
}

void Base_Variable::assign_distr_scalar_arg(const std::string &arg_name, double x)
{
    _distr_ptr_->assign_scalar_arg(arg_name, x);
}
void Base_Variable::assign_distr_pointer_arg(const std::string &arg_name, const double *y)
{
    _distr_ptr_->assign_pointer_arg(arg_name, y);
}

/* Assign reshaping to distribution. */
void Base_Variable::assign_distr_range(double lower, double upper)
{
    _distr_ptr_->assign_reshaping(_is_discrete_, lower, upper);
}

/* Assign proposer. */
void Base_Variable::assign_proposer(const std::string &rule_name,
                                    double rule_scale,
                                    double rule_from)
{
    if (tools::is_na(rule_name))
    {
        // if no rule is defined constant variable
        _is_constant_ = true;
        _func_propose_ = &Base_Variable::_propose_error_;
        _prop_ptr_ = proposer_switchboard();
    }
    else
    {
        // variable to sample
        _is_constant_ = false;

        // define proposer rule with functional programming
        if (rule_name == "gibbs")
        {
            //----- Gibbs proposer rule -----//
            // assign the pointer to the proposer
            _prop_ptr_ = std::make_shared<Gibbs_Proposer>(_distr_ptr_);
            // sampler
            _func_propose_ = &Base_Variable::_propose_gibbs_sample_;
            // logratio
            _func_move_logratio_ = &Base_Variable::_move_logratio_gibbs_;
            // validation
            _func_valid_move_ = &Base_Variable::_valid_move_gibbs_;
        }
        else
        {
            //----- Other proposer rules -----//
            // value from which to start proposal
            if ((!tools::is_na(rule_from)) && (!_distr_ptr_->valid(rule_from)))
            {
                std::ostringstream oss;
                oss << "Proposer rule from value " << rule_from << " out of range"
                    << " {min:" << tools::natod(show_range_min())
                    << ", max:" << tools::natod(show_range_max()) << "}."
                    << " If intentional, ignore this message.";
                tools::warning(oss.str());
            }

            // assign the pointer to the proposer
            _prop_ptr_ = proposer_switchboard(rule_name, rule_scale, rule_from, _is_discrete_);
            // sampler
            if (is_from_zero())
                _func_propose_ = &Base_Variable::_propose_from_zero_;
            else
                _func_propose_ = &Base_Variable::_propose_from_previous_;
        }
    }
}

/* Store the buffer and change value. */
void Base_Variable::set_value(double x)
{
    buffer = value;
    value = x;
}

/* Restore buffered value. */
void Base_Variable::restore_value() { value = buffer; }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Base_Variable - METROPOLIS-HASTINGS  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Variable::get_move_logratio() const { return (this->*_func_move_logratio_)(); }
bool Base_Variable::valid_move() const { return (this->*_func_valid_move_)(); }

//----- distr_ptr : prior/model methods -----//
void Base_Variable::rand(rnd_gen &prng) { set_value(_distr_ptr_->rand(prng)); }

double Base_Variable::distr_logpdf() const { return _distr_ptr_->logpdf(value); }
double Base_Variable::distr_logratio() const { return _distr_ptr_->logratio(buffer, value); }

//----- prop_ptr : propose methods -----//
void Base_Variable::propose(rnd_gen &prng) { return set_value((this->*_func_propose_)(prng)); }
double Base_Variable::propose_logratio() const { return (_prop_ptr_->logratio)(buffer, value); }

//----- prop_ptr : adaptive methods -----//
double Base_Variable::get_accept_rate() const { return _prop_ptr_->get_accept_rate(); }
void Base_Variable::update_accept_rate(bool is_accepted) { _prop_ptr_->update_accept_rate(is_accepted); }
void Base_Variable::reset_accept_rate() { _prop_ptr_->reset_accept_rate(); }
void Base_Variable::adapt_rule(double ratio) { (_prop_ptr_->adapt)(ratio); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base_Variable - SHOW  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Variable::get_distr_name() const { return _distr_ptr_->get_name(); }
const std::string Base_Variable::show_distr_arg_names() const { return _distr_ptr_->show_arg_names(); }
const std::string Base_Variable::show_distr_arg_values() const { return _distr_ptr_->show_arg_values(); }

const std::string &Base_Variable::get_prop_name() const { return _prop_ptr_->get_name(); }
const std::string Base_Variable::show_prop_arg_names() const { return _prop_ptr_->show_arg_names(); }
const std::string Base_Variable::show_prop_arg_values() const { return _prop_ptr_->show_arg_values(); }
const std::string Base_Variable::show_prop_scale() const { return _prop_ptr_->show_scale(); }

const std::string Base_Variable::show_range_min() const { return _distr_ptr_->show_range_min(); }
const std::string Base_Variable::show_range_max() const { return _distr_ptr_->show_range_max(); }

const std::string Base_Variable::show() const
{
    std::ostringstream oss;
    oss << name
        << C_DELIM << tools::natod(get_value())
        << C_DELIM << tools::truth(_is_constant_)
        << C_DELIM << tools::truth(_is_discrete_)
        << C_DELIM << get_distr_name()
        << C_DELIM << show_distr_arg_names()
        << C_DELIM << show_distr_arg_values()
        << C_DELIM << get_prop_name()
        << C_DELIM << tools::natod(show_prop_scale())
        << C_DELIM << tools::natod(show_range_min())
        << C_DELIM << tools::natod(show_range_max());

    return oss.str();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base_Variable - functional  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== Valid move : functional methods =====//
bool Base_Variable::_valid_move_default_() const { return _distr_ptr_->valid(value); }
bool Base_Variable::_valid_move_gibbs_() const { return true; }

//===== Move logratio : functional methods =====//
double Base_Variable::_move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
double Base_Variable::_move_logratio_gibbs_() const { return 0.0; }

//===== Propose : functional methods =====//
double Base_Variable::_propose_error_(rnd_gen &prng)
{
    UNUSED(prng);
    throw std::invalid_argument(
        "Call to proposer of constant variable '" + get_name() + "'.");
    return 0.;
}
double Base_Variable::_propose_gibbs_sample_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, 0); }
double Base_Variable::_propose_from_previous_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, get_value()); }
double Base_Variable::_propose_from_zero_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, get_from_value()); }
