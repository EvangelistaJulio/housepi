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
    : pDistr_(distr)
{
    // default proposer name derives from distribution
    proposer_name_ = pDistr_->name();

    // set to discrete
    pDistr_->init_reshaping(is_discrete);

    // assign user rule scale
    if (!tools::is_na(rule_scale))
        user_rule_scale_ = std::fabs(rule_scale);

    // by default proposal is from previous state
    is_start_from_zero_ = false;

    // proposal alwats starting by given value
    if (!tools::is_na(rule_start_value))
    {
        is_start_from_zero_ = true;
        set_rule_start_value(rule_start_value);
    }

    // load which are the parameter to be adapted
    n_args_to_adapt_ = args_to_adapt.size();
    args_to_adapt_list_ = new ushort[n_args_to_adapt_];

    for (ushort c = 0; c < n_args_to_adapt_; ++c)
        args_to_adapt_list_[c] = args_to_adapt.at(c);
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
    for (ushort c = 0; c < n_args_to_adapt_; ++c)
    {
        ushort arg_id = args_to_adapt_list_[c];
        (pDistr_->set_arg)(arg_id, factor * (pDistr_->arg)(arg_id));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Base Show - adaptation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

const std::string Proposer::str_scale() const
{
    if (n_args_to_adapt_ == 0)
        return HSP_NA_STRING;
    else
        return std::to_string(std::fabs((pDistr_->arg)(args_to_adapt_list_[0])));
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
        pDistr_ = distr;

        proposer_name_ = "draw";
        // proposer functional method
        propose_ = &Draw_Proposer::do_propose_draw_rand_;
        // logratio functional method
        logratio_ = &Draw_Proposer::do_logratio_pdf_;
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
    Uniform_Proposer(double rule_scale = HSP_NA_DOUBLE,
                     double rule_from = HSP_NA_DOUBLE,
                     bool is_discrete = false)
        : Proposer("unif", rule_scale, rule_from, {0, 1}, is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        propose_ = &Uniform_Proposer::do_propose_sum_to_previous_;
        // logratio functional method
        logratio_ = &Uniform_Proposer::do_logratio_symm_;
    }
    ~Uniform_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        pDistr_->init_attr("min", HSP_ATTR_ARGUMENT, -user_rule_scale_);
        pDistr_->init_attr("max", HSP_ATTR_ARGUMENT, user_rule_scale_);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Normal_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Proposer
{
public:
    //===== CONSTRUCTOR =====//
    Normal_Proposer(double rule_scale = HSP_NA_DOUBLE,
                    double rule_from = HSP_NA_DOUBLE,
                    bool is_discrete = false)
        : Proposer("norm", rule_scale, rule_from, {1}, is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        propose_ = &Normal_Proposer::do_propose_sum_to_previous_;
        // logratio functional method
        if (is_rule_start_from_zero())
            logratio_ = &Normal_Proposer::do_logratio_sum_to_zero_;
        else
            logratio_ = &Normal_Proposer::do_logratio_symm_;
    }
    ~Normal_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        pDistr_->init_attr("mean", HSP_ATTR_ARGUMENT, 0.);
        pDistr_->init_attr("sd", HSP_ATTR_ARGUMENT, user_rule_scale_);
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
                     double rule_from = HSP_NA_DOUBLE,
                     bool is_discrete = false)
        : Proposer("lnorm", rule_scale, rule_from, {1}, is_discrete)
    {
        // initialize
        user_reset();
        // proposer functional method
        propose_ = &Lognorm_Proposer::do_propose_times_previous_;
        // logratio functional method
        if (is_rule_start_from_zero())
            logratio_ = &Lognorm_Proposer::do_logratio_times_zero_;
        else
            logratio_ = &Lognorm_Proposer::do_logratio_lognorm_symm_;
    }
    ~Lognorm_Proposer() = default;

    //===== SETTER =====//
    void user_reset() override final
    {
        pDistr_->init_attr("meanlog", HSP_ATTR_ARGUMENT, 0.);
        pDistr_->init_attr("sdlog", HSP_ATTR_ARGUMENT, user_rule_scale_);
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Switchboard to select the custom metropolis proposer rule. */
std::shared_ptr<Proposer> proposer_switch(const std::string &rule_name = HSP_NA_STRING,
                                          double rule_scale = HSP_NA_DOUBLE,
                                          double rule_from = HSP_NA_DOUBLE,
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
        return_value_ = &Variable::do_return_error_;
        get_ = &Variable::do_get_default_;
        get_ptr_ = &Variable::do_get_ptr_default_;

        variable_name_ = name;
        user_value_ = value;
        value_ = user_value_;
        is_discrete_ = is_discrete;
        is_constant_ = is_constant;
    }
    else
    {
        return_value_ = &Variable::do_return_func_;
        get_ = &Variable::do_get_func_return_;
        get_ptr_ = &Variable::do_get_ptr_func_return_;
        pFunc_ = compose_func_registries(
            func_name, "function",
            {&UserFuncRegistry::inst(),
             &TranscendentalRegistry::inst(),
             &ArithmeticRegistry::inst(),
             &ContinuosDistrRegistry::inst(),
             &DiscreteDistrRegistry::inst(),
             &InfecProfileRegistry::inst()});
        is_constant_ = false;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Variable - Setter (Distribution)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign distribution with constant arguments. */
void Variable::init_distr(const std::string &distr_name)
{
    pDistr_ = compose_distr_registries(
        distr_name,
        "distribution",
        {&ContinuosDistrRegistry::inst(),
         &DiscreteDistrRegistry::inst()});
}
/* Initialize distribution reshaping. */
void Variable::init_distr_reshaping()
{
    pDistr_->init_reshaping(is_discrete_);
    // check that non-constant variable has non-zero dimensional range
    if ((!is_constant_) && (pDistr_->get_min() == pDistr_->get_max()))
        throw std::invalid_argument(
            "Degenerate range " + pDistr_->str_range() + " defined for non constant variable \"" + name() + "\".");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Variable - Setter (Proposer)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign proposer. */
void Variable::init_proposer(const std::string &rule_name,
                             double rule_scale,
                             double rule_from)
{
    if (tools::is_na(rule_name))
    {
        // if no rule is defined constant variable
        is_constant_ = true;
        propose_ = &Variable::do_propose_error_;
        pProp_ = proposer_switch();
    }
    else
    {
        // variable to sample
        is_constant_ = false;

        // define proposer rule with functional programming
        if (rule_name == "draw")
        {
            //----- Draw proposer rule -----//
            // assign the pointer to the proposer
            pProp_ = std::make_shared<Draw_Proposer>(pDistr_);

            // sampler
            propose_ = &Variable::do_propose_draw_sample_;
            // logratio
            move_logratio_ = &Variable::do_move_logratio_draw_;
            // validation
            valid_move_ = &Variable::do_valid_move_draw_;
        }
        else
        {
            //----- Other proposer rules -----//
            // value from which to start proposal
            if ((!tools::is_na(rule_from)) && (!pDistr_->valid(rule_from)))
                tools::warning(
                    "Proposer rule from value " + std::to_string(rule_from) + " out of range" + " {min:" + str_min() + ", max:" + str_max() + "}.");

            // assign the pointer to the proposer
            pProp_ = proposer_switch(rule_name, rule_scale, rule_from, is_discrete());

            // raise a warning if "lnorm" rule when varaible can be negative
            if ((pDistr_->get_min() < 0) && (rule_name == "lnorm"))
                tools::warning(
                    "Lognormal proposer rule for non-positive variable \"" + variable_name_ + "\" .");

            // sampler
            if (is_from_zero())
                propose_ = &Variable::do_propose_from_zero_;
            else
                propose_ = &Variable::do_propose_from_previous_;
        }
    }
}
