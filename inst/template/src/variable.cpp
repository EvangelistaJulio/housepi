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
///////////////////////////////////  Proposer : constructor  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Proposer constructor from Distribution pointer. */
Proposer::Proposer(std::shared_ptr<Distribution> distr,
                   double rule_scale,
                   double rule_start_value,
                   const std::vector<ushort> &args_to_adapt,
                   bool is_discrete)
    : _distr_ptr(distr)
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

/* Proposer constructor from Distribution name. */
Proposer::Proposer(const std::string &name,
                   double rule_scale,
                   double rule_start_value,
                   const std::vector<ushort> &args_to_adapt,
                   bool is_discrete)
    : Proposer(ContinuosDistrRegistry::inst().get_distr(name),
               rule_scale, rule_start_value, args_to_adapt, is_discrete) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Base Proposer - adaptation  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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

const std::string Proposer::show_scale() const
{
    if (_n_args_to_adapt == 0)
        return HOUSEPI_NA_STRING;
    else
        return std::to_string(std::fabs((_distr_ptr->get_arg)(_args_to_adapt_list[0])));
}

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
        _do_propose = &Draw_Proposer::__propose_draw_rand;
        // logratio functional method
        _do_logratio = &Draw_Proposer::__logratio_pdf;
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
    Uniform_Proposer(double rule_scale = HOUSEPI_NA_DOUBLE,
                     double rule_from = HOUSEPI_NA_DOUBLE,
                     bool is_discrete = false)
        : Proposer("unif", rule_scale, rule_from, {0, 1}, is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        _do_propose = &Uniform_Proposer::__propose_sum_to_previous;
        // logratio functional method
        _do_logratio = &Uniform_Proposer::__logratio_symm;
    }
    ~Uniform_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr->assign_attr("min", ATTR_ARGUMENT, -_user_rule_scale);
        _distr_ptr->assign_attr("max", ATTR_ARGUMENT, _user_rule_scale);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Normal_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Normal_Proposer(double rule_scale = HOUSEPI_NA_DOUBLE,
                    double rule_from = HOUSEPI_NA_DOUBLE,
                    bool is_discrete = false)
        : Proposer("norm", rule_scale, rule_from, {1}, is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        _do_propose = &Normal_Proposer::__propose_sum_to_previous;
        // logratio functional method
        if (is_rule_start_from_zero())
            _do_logratio = &Normal_Proposer::__logratio_sum_to_zero;
        else
            _do_logratio = &Normal_Proposer::__logratio_symm;
    }
    ~Normal_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr->assign_attr("mean", ATTR_ARGUMENT, 0.);
        _distr_ptr->assign_attr("sd", ATTR_ARGUMENT, _user_rule_scale);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Lognorm_Proposer(double rule_scale = 1.,
                     double rule_from = HOUSEPI_NA_DOUBLE,
                     bool is_discrete = false)
        : Proposer("lnorm", rule_scale, rule_from, {1}, is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        _do_propose = &Lognorm_Proposer::__propose_times_previous;
        // logratio functional method
        if (is_rule_start_from_zero())
            _do_logratio = &Lognorm_Proposer::__logratio_times_zero;
        else
            _do_logratio = &Lognorm_Proposer::__logratio_lognorm_symm;
    }
    ~Lognorm_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        _distr_ptr->assign_attr("meanlog", ATTR_ARGUMENT, 0.);
        _distr_ptr->assign_attr("sdlog", ATTR_ARGUMENT, _user_rule_scale);
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Switchboard to select the custom metropolis proposer rule. */
std::shared_ptr<Proposer> proposer_switch(const std::string &rule_name = HOUSEPI_NA_STRING,
                                          double rule_scale = HOUSEPI_NA_DOUBLE,
                                          double rule_from = HOUSEPI_NA_DOUBLE,
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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable - Setter  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Restore user defined variable initialization. */
void Variable::init(const std::string &name, double value,
                    bool is_discrete, bool is_constant, const std::string &func_name)
{
    if (tools::is_na(func_name))
    {
        _do_get = &Variable::__get_default;
        _do_get_pointer = &Variable::__get_pointer_default;

        _variable_name = name;
        _user_value = value;
        _value = _user_value;
        _is_discrete = is_discrete;
        _is_constant = is_constant;
    }
    else
    {
        _do_get = &Variable::__get_func_return;
        _do_get_pointer = &Variable::__get_pointer_func_return;
        _func_ptr = compose_func_registries(
            func_name,
            "function",
            {&ArithmeticRegistry::inst(),
             &TranscendentalRegistry::inst(),
             &ContinuosDistrRegistry::inst(),
             &DiscreteDistrRegistry::inst(),
             &InfecProfileRegistry::inst()});
        _is_constant = false;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Variable - Setter (Distribution)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign distribution with constant arguments. */
void Variable::assign_distr(const std::string &distr_name)
{
    _distr_ptr = compose_distr_registries(
        distr_name,
        "distribution",
        {&ContinuosDistrRegistry::inst(),
         &DiscreteDistrRegistry::inst()});
}
/* Initialize distribution reshaping. */
void Variable::init_distr_reshaping()
{
    _distr_ptr->init_reshaping(_is_discrete);
    // check that non-constant variable has non-zero dimensional range
    if ((!_is_constant) && (_distr_ptr->get_min() == _distr_ptr->get_max()))
        throw std::invalid_argument(
            "Degenerate range " + _distr_ptr->show_range() + " defined for non constant variable \"" + show_name() + "\".");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Variable - Setter (Proposer)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign proposer. */
void Variable::assign_proposer(const std::string &rule_name,
                               double rule_scale,
                               double rule_from)
{
    if (tools::is_na(rule_name))
    {
        // if no rule is defined constant variable
        _is_constant = true;
        _do_propose = &Variable::__propose_error;
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
            _do_propose = &Variable::__propose_draw_sample;
            // logratio
            _do_move_logratio = &Variable::__move_logratio_draw;
            // validation
            _do_valid_move = &Variable::__valid_move_draw;
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
                _do_propose = &Variable::__propose_from_zero;
            else
                _do_propose = &Variable::__propose_from_previous;
        }
    }
}
