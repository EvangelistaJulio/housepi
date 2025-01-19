/*
 * variable.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "variable.h"
#include "distribution.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Proposer : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Base class for Metropolis proposers. */
class Base_Proposer
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//
    Base_Proposer() {};
    Base_Proposer(std::shared_ptr<Base_Distribution>,
                  double rule_scale = 1.0,
                  double rule_start_value = NULL_DOUBLE,
                  const std::vector<ushort> &args_to_adapt = {},
                  bool is_discrete = false);
    virtual ~Base_Proposer() {};

    //==================//
    //===== GETTER =====//
    //==================//

    const std::string &get_name() const { return this->name; }
    double get_arg(ushort argID) const { return (_distr_ptr_->get_arg)(argID); }
    double get_from_value() const { return from_value; }
    bool is_rule_start_from_zero() const { return _is_start_from_zero_; }
    double propose(rnd_gen &prng, double start = 0.0) { return (this->*_func_propose_)(prng, start); }
    double logratio(double xold, double xnew) const { return (this->*_func_logratio_)(xold, xnew); }
    double rand(rnd_gen &prng) { return (_distr_ptr_->rand)(prng); }

    //==================//
    //===== SETTER =====//
    //==================//

    virtual void user_reset() {};
    void set_rule_start_value(double x) { from_value = x; }

    //======================//
    //===== ADAPTATION =====//
    //======================//

    double get_accept_rate() const { return accept_rate; }
    void update_accept_rate(bool is_accepted)
    {
        accept_rate = (double)(accept_rate * n_proposed_moves + is_accepted) / (n_proposed_moves + 1);
        n_proposed_moves += 1;
    }
    void reset_accept_rate()
    {
        accept_rate = 0.0;
        n_proposed_moves = 0;
    };

    // increase spread if acceptation rate is high and viceversa
    void adapt(double ratio)
    {
        double factor = 1. + (1. - ratio) * (get_accept_rate() - TARGET_ACCEPT_RATE);
        // loop over all arguments to adapt
        for (const ushort &argID : _args_to_adapt_)
            (_distr_ptr_->set_arg)(argID,
                                   (_distr_ptr_->get_arg)(argID)*factor);
    }

    //================//
    //===== SHOW =====//
    //================//

    const std::string show_arg_names() const { return _distr_ptr_->show_arg_names(); }
    const std::string show_arg_values() const { return _distr_ptr_->show_arg_values(); }
    const std::string show_scale() const
    {
        if (_args_to_adapt_.empty())
            return NULL_STRING;
        else
            return std::to_string(std::fabs((_distr_ptr_->get_arg)(_args_to_adapt_[0])));
    }

protected:
    std::string name = NULL_STRING;
    std::shared_ptr<Base_Distribution> _distr_ptr_;

    double _user_rule_scale_ = 1.0;    // scale of the Metropolis rule
    bool _is_start_from_zero_ = false; // wheter to move from given value
    double from_value = NULL_DOUBLE;   // given value where to start from
    double accept_rate = 0.0;          // acceptance rate
    ulong n_proposed_moves = 0;        // number of proposed moves

    std::vector<ushort> _args_to_adapt_ = {}; // which arguments should be adapted

    //===== propose : functional methods =====//
    double (Base_Proposer::*_func_propose_)(rnd_gen &, double) = &Base_Proposer::_propose_undefined_;
    double _propose_undefined_(rnd_gen &prng, double start);
    double _propose_sum_to_previous_(rnd_gen &prng, double start);
    double _propose_times_previous_(rnd_gen &prng, double start);
    double _propose_draw_rand_(rnd_gen &prng, double start);

    //===== logratio : functional methods =====//
    double (Base_Proposer::*_func_logratio_)(double, double) const = &Base_Proposer::_logratio_undefined_;
    double _logratio_undefined_(double xold, double xnew) const;
    double _logratio_pdf_(double xold, double xnew) const;
    double _logratio_sum_to_zero_(double xold, double xnew) const;
    double _logratio_times_zero_(double xold, double xnew) const;
    double _logratio_symm_(double xold, double xnew) const;
    double _logratio_lognorm_symm_(double xold, double xnew) const;
};

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
    _distr_ptr_->init_reshaping(is_discrete);

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
    Gibbs_Proposer(std::shared_ptr<Base_Distribution> distr)
    {
        _distr_ptr_ = distr;

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

/* Switchboard to select the custom metropolis proposer rule. */
std::shared_ptr<Base_Proposer> proposer_switchboard(const std::string &rule_name = NULL_STRING,
                                                    double rule_scale = NULL_DOUBLE,
                                                    double rule_from = NULL_DOUBLE,
                                                    bool is_discrete = false)
{
    const str_vec options = {"gibbs", "unif", "norm", "lnorm"};

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
                "Unknown proposer rule \"" + rule_name + "\"." + tools::find_closest(rule_name, options));
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

class Base_Variable::Impl
{
public:
    std::shared_ptr<Base_Distribution> _distr_ptr_; // prior/model distribution of the variable
    std::shared_ptr<Base_Proposer> _prop_ptr_;      // proposer for Monte Carlo

    std::string name = NULL_STRING; // name of the variable
    double value = NULL_DOUBLE;     // value of the variable
    double buffer = NULL_DOUBLE;    // buffer of the variable

    double _user_value_ = NULL_DOUBLE; // initial user defined value
    bool _is_discrete_ = false;        // wheter discrete or not
    bool _is_constant_ = true;         // wheter constant or not

    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Impl() : _distr_ptr_(std::make_shared<Base_Distribution>()),
             _prop_ptr_(std::make_shared<Base_Proposer>()) {}
    ~Impl() = default;

    //===============================//
    //===== METROPOLIS-HASTINGS =====//
    //===============================//

    void set_value(double x);
    void restore_value();

    double get_move_logratio() const;
    bool valid_move() const;

    //----- distr_ptr : prior/model methods -----//
    void rand(rnd_gen &prng);

    double distr_logpdf() const;
    double distr_logratio() const;

    //----- prop_ptr : propose methods -----//
    void propose(rnd_gen &prng);
    double propose_logratio() const;

    //----- prop_ptr : adaptive methods -----//
    double get_accept_rate() const;
    void update_accept_rate(bool is_accepted);
    void reset_accept_rate();
    void adapt_rule(double ratio);

    //===== Valid move : functional methods =====//
    bool (Impl::*_func_valid_move_)() const = &Impl::_valid_move_default_;
    bool _valid_move_default_() const;
    bool _valid_move_gibbs_() const;

    //===== Move logratio : functional methods =====//
    double (Impl::*_func_move_logratio_)() const = &Impl::_move_logratio_default_;
    double _move_logratio_default_() const;
    double _move_logratio_gibbs_() const;

    //===== Propose : functional methods =====//
    double (Impl::*_func_propose_)(rnd_gen &) = &Impl::_propose_error_;
    double _propose_error_(rnd_gen &prng);
    double _propose_gibbs_sample_(rnd_gen &prng);
    double _propose_from_previous_(rnd_gen &prng);
    double _propose_from_zero_(rnd_gen &prng);
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Base_Variable::Impl - functional  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== Valid move : functional methods =====//
bool Base_Variable::Impl::_valid_move_default_() const { return _distr_ptr_->valid(value); }
bool Base_Variable::Impl::_valid_move_gibbs_() const { return true; }

//===== Move logratio : functional methods =====//
double Base_Variable::Impl::_move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
double Base_Variable::Impl::_move_logratio_gibbs_() const { return 0.0; }

//===== Propose : functional methods =====//
double Base_Variable::Impl::_propose_error_(rnd_gen &prng)
{
    UNUSED(prng);
    throw std::invalid_argument(
        "Call to proposer of constant variable \"" + name + "\".");
    return 0.;
}
double Base_Variable::Impl::_propose_gibbs_sample_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, 0); }
double Base_Variable::Impl::_propose_from_previous_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, value); }
double Base_Variable::Impl::_propose_from_zero_(rnd_gen &prng) { return (_prop_ptr_->propose)(prng, _prop_ptr_->get_from_value()); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base_Variable::Impl - METROPOLIS-HASTINGS  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store the buffer and change value. */
inline void Base_Variable::Impl::set_value(double x)
{
    buffer = value;
    value = x;
}

/* Restore buffered value. */
inline void Base_Variable::Impl::restore_value() { value = buffer; }

inline double Base_Variable::Impl::get_move_logratio() const { return (this->*_func_move_logratio_)(); }
inline bool Base_Variable::Impl::valid_move() const { return (this->*_func_valid_move_)(); }

//----- distr_ptr : prior/model methods -----//
inline void Base_Variable::Impl::rand(rnd_gen &prng) { set_value(_distr_ptr_->rand(prng)); }

inline double Base_Variable::Impl::distr_logpdf() const { return _distr_ptr_->logpdf(value); }
inline double Base_Variable::Impl::distr_logratio() const { return _distr_ptr_->logratio(buffer, value); }

//----- prop_ptr : propose methods -----//
inline void Base_Variable::Impl::propose(rnd_gen &prng) { set_value((this->*_func_propose_)(prng)); }

inline double Base_Variable::Impl::propose_logratio() const { return (_prop_ptr_->logratio)(buffer, value); }

//----- prop_ptr : adaptive methods -----//
inline double Base_Variable::Impl::get_accept_rate() const { return _prop_ptr_->get_accept_rate(); }
inline void Base_Variable::Impl::update_accept_rate(bool is_accepted) { _prop_ptr_->update_accept_rate(is_accepted); }
inline void Base_Variable::Impl::reset_accept_rate() { _prop_ptr_->reset_accept_rate(); }
inline void Base_Variable::Impl::adapt_rule(double ratio) { (_prop_ptr_->adapt)(ratio); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base_Variable - CONSTRUCTOR  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Base_Variable::Base_Variable(const std::string &name,
                             double value,
                             bool is_discrete) : pImpl(std::make_unique<Impl>())
{
    init(name, value, is_discrete);
}
Base_Variable::Base_Variable(Base_Variable &&other) noexcept = default;
Base_Variable &Base_Variable::operator=(Base_Variable &&other) noexcept = default;
Base_Variable::~Base_Variable() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Variable - GETTER  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Variable::get_name() const { return pImpl->name; }
double Base_Variable::get_value() const { return pImpl->value; }
const double *Base_Variable::get_pointer() const { return &pImpl->value; }
bool Base_Variable::is_constant() const { return pImpl->_is_constant_; }
bool Base_Variable::is_discrete() const { return pImpl->_is_discrete_; }
bool Base_Variable::is_from_zero() const { return (pImpl->_prop_ptr_->is_rule_start_from_zero)(); }
double Base_Variable::get_from_value() const { return (pImpl->_prop_ptr_->get_from_value)(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Variable - SETTER  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Restore user defined variable initialization. */
void Base_Variable::init(const std::string &name, double value, bool is_discrete)
{
    pImpl->name = name;
    pImpl->_user_value_ = value;
    pImpl->value = value;
    pImpl->_is_discrete_ = is_discrete;
}

/* Restore user defined variable initialization. */
void Base_Variable::user_reset()
{
    pImpl->value = pImpl->_user_value_;
    pImpl->_prop_ptr_->user_reset();
}

/* Assign distribution with constant arguments. */
void Base_Variable::assign_distribution(const std::string &distr_name)
{
    pImpl->_distr_ptr_ = distribution_switchboard(distr_name);
}

void Base_Variable::assign_distr_scalar_arg(const std::string &arg_name, double x)
{
    pImpl->_distr_ptr_->assign_scalar_arg(arg_name, x);
}
void Base_Variable::assign_distr_pointer_arg(const std::string &arg_name, const double *y)
{
    pImpl->_distr_ptr_->assign_pointer_arg(arg_name, y);
}

void Base_Variable::assign_distr_scalar_range_min(double x)
{
    pImpl->_distr_ptr_->assign_scalar_range_min(x);
}
void Base_Variable::assign_distr_pointer_range_min(const double *y)
{
    pImpl->_distr_ptr_->assign_pointer_range_min(y);
}

void Base_Variable::assign_distr_scalar_range_max(double x)
{
    pImpl->_distr_ptr_->assign_scalar_range_max(x);
}
void Base_Variable::assign_distr_pointer_range_max(const double *y)
{
    pImpl->_distr_ptr_->assign_pointer_range_max(y);
}

void Base_Variable::init_distr_reshaping()
{
    pImpl->_distr_ptr_->init_reshaping(is_discrete());
}

/* Assign proposer. */
void Base_Variable::assign_proposer(const std::string &rule_name,
                                    double rule_scale,
                                    double rule_from)
{
    if (tools::is_na(rule_name))
    {
        // if no rule is defined constant variable
        pImpl->_is_constant_ = true;
        pImpl->_func_propose_ = &Impl::_propose_error_;
        pImpl->_prop_ptr_ = proposer_switchboard();
    }
    else
    {
        // variable to sample
        pImpl->_is_constant_ = false;

        // define proposer rule with functional programming
        if (rule_name == "gibbs")
        {
            //----- Gibbs proposer rule -----//
            // assign the pointer to the proposer
            pImpl->_prop_ptr_ = std::make_shared<Gibbs_Proposer>(pImpl->_distr_ptr_);

            // sampler
            pImpl->_func_propose_ = &Impl::_propose_gibbs_sample_;
            // logratio
            pImpl->_func_move_logratio_ = &Impl::_move_logratio_gibbs_;
            // validation
            pImpl->_func_valid_move_ = &Impl::_valid_move_gibbs_;
        }
        else
        {
            //----- Other proposer rules -----//
            // value from which to start proposal
            if ((!tools::is_na(rule_from)) && (!pImpl->_distr_ptr_->valid(rule_from)))
            {
                std::ostringstream oss;
                oss << "Proposer rule from value " << rule_from << " out of range"
                    << " {min:" << tools::natod(show_range_min())
                    << ", max:" << tools::natod(show_range_max()) << "}."
                    << " If intentional, ignore this message.";
                tools::warning(oss.str());
            }

            // assign the pointer to the proposer
            pImpl->_prop_ptr_ = proposer_switchboard(rule_name, rule_scale, rule_from, is_discrete());
            // sampler
            if (is_from_zero())
                pImpl->_func_propose_ = &Impl::_propose_from_zero_;
            else
                pImpl->_func_propose_ = &Impl::_propose_from_previous_;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Base_Variable - METROPOLIS-HASTINGS  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store the buffer and change value. */
void Base_Variable::set_value(double x) { pImpl->set_value(x); }

/* Restore buffered value. */
void Base_Variable::restore_value() { pImpl->restore_value(); }

double Base_Variable::get_move_logratio() const { return pImpl->get_move_logratio(); }
bool Base_Variable::valid_move() const { return pImpl->valid_move(); }

//----- distr_ptr : prior/model methods -----//
void Base_Variable::rand(rnd_gen &prng) { pImpl->rand(prng); }

double Base_Variable::distr_logpdf() const { return pImpl->distr_logpdf(); }
double Base_Variable::distr_logratio() const { return pImpl->distr_logratio(); }

//----- prop_ptr : propose methods -----//
void Base_Variable::propose(rnd_gen &prng) { pImpl->propose(prng); }
double Base_Variable::propose_logratio() const { return pImpl->propose_logratio(); }

//----- prop_ptr : adaptive methods -----//
double Base_Variable::get_accept_rate() const { return pImpl->get_accept_rate(); }
void Base_Variable::update_accept_rate(bool is_accepted) { pImpl->update_accept_rate(is_accepted); }
void Base_Variable::reset_accept_rate() { pImpl->reset_accept_rate(); }
void Base_Variable::adapt_rule(double ratio) { pImpl->adapt_rule(ratio); }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Base_Variable - SHOW  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string &Base_Variable::get_distr_name() const { return pImpl->_distr_ptr_->get_name(); }
const std::string Base_Variable::show_distr_arg_names() const { return pImpl->_distr_ptr_->show_arg_names(); }
const std::string Base_Variable::show_distr_arg_values() const { return pImpl->_distr_ptr_->show_arg_values(); }

const std::string &Base_Variable::get_prop_name() const { return pImpl->_prop_ptr_->get_name(); }
const std::string Base_Variable::show_prop_arg_names() const { return pImpl->_prop_ptr_->show_arg_names(); }
const std::string Base_Variable::show_prop_arg_values() const { return pImpl->_prop_ptr_->show_arg_values(); }
const std::string Base_Variable::show_prop_scale() const { return pImpl->_prop_ptr_->show_scale(); }

const std::string Base_Variable::show_range_min() const { return pImpl->_distr_ptr_->show_range_min(); }
const std::string Base_Variable::show_range_max() const { return pImpl->_distr_ptr_->show_range_max(); }

const std::string Base_Variable::show() const
{
    std::ostringstream oss;
    oss << pImpl->name
        << C_DELIM << tools::natod(get_value())
        << C_DELIM << tools::truth(is_constant())
        << C_DELIM << tools::truth(is_discrete())
        << C_DELIM << get_distr_name()
        << C_DELIM << show_distr_arg_names()
        << C_DELIM << show_distr_arg_values()
        << C_DELIM << get_prop_name()
        << C_DELIM << tools::natod(show_prop_scale())
        << C_DELIM << tools::natod(show_range_min())
        << C_DELIM << tools::natod(show_range_max());

    return oss.str();
}
