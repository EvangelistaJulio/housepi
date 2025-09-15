/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * variable.cpp
 *
 */

#include "variable.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Proposer : constructor  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Proposer constructor from Distribution pointer.
Proposer::Proposer(std::shared_ptr<Distribution> distr,
                   double rule_scale,
                   double rule_start_value,
                   const std::vector<unsigned short> &args_to_adapt,
                   bool is_discrete)
    : pDistr_(distr)
{
    // default proposer name derives from distribution
    name_ = pDistr_->name();

    // set to discrete
    if (is_discrete)
        pDistr_->init_discretization();
    // set range and truncation
    pDistr_->init_truncation();

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
    args_to_adapt_list_ = new unsigned short[n_args_to_adapt_];

    for (unsigned short c = 0; c < n_args_to_adapt_; ++c)
        args_to_adapt_list_[c] = args_to_adapt[c];
};

/// Proposer constructor from Distribution name.
Proposer::Proposer(const std::string &name,
                   double rule_scale,
                   double rule_start_value,
                   const std::vector<unsigned short> &args_to_adapt,
                   bool is_discrete)
    : Proposer(ContinuosDistrRegistry::inst().get_distr(name),
               rule_scale, rule_start_value, args_to_adapt, is_discrete) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Proposer - adaptation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Increase spread if acceptation rate is high and viceversa
void Proposer::adapt(double ratio)
{
    const double factor = 1. + (1. - ratio) * (get_accept_rate() - HSP_TARGET_ACCEPT_RATE);
    // loop over all arguments to adapt
    for (unsigned short c = 0; c < n_args_to_adapt_; ++c)
    {
        unsigned short arg_id = args_to_adapt_list_[c];
        (pDistr_->set_arg)(arg_id, factor * (pDistr_->arg)(arg_id));
    }
}

/// Format proposer scale as a string.
std::string Proposer::str_scale() const
{
    if (n_args_to_adapt_ == 0)
        return HSP_NA_STRING;
    else
        return std::to_string(std::fabs((pDistr_->arg)(args_to_adapt_list_[0])));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// *** Sample_Proposer : class *** /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Sample_Proposer : public Proposer
{
public:
    Sample_Proposer(std::shared_ptr<Distribution> distr)
    {
        name_ = "sample";
        propose_ = &Sample_Proposer::do_propose_sample_rand_; // proposer functional method
        logratio_ = &Sample_Proposer::do_logratio_pdf_;       // logratio functional method
        pDistr_ = distr;                                      // link to model distribution
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Uniform_Proposer : class *** /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Proposer
{
public:
    Uniform_Proposer(double rule_scale = HSP_NA_DOUBLE,
                     double rule_from = HSP_NA_DOUBLE,
                     bool is_discrete = false)
        : Proposer("unif", rule_scale, rule_from, {0, 1}, is_discrete)
    {
        user_reset();                                              // initialize
        propose_ = &Uniform_Proposer::do_propose_sum_to_previous_; // proposer functional method
        logratio_ = &Uniform_Proposer::do_logratio_symm_;          // logratio functional method
    }

    void user_reset() override final
    {
        pDistr_->init_attr("min", HSP_ATTR_ARGUMENT, -user_rule_scale_);
        pDistr_->init_attr("max", HSP_ATTR_ARGUMENT, user_rule_scale_);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Normal_Proposer : class *** //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Proposer
{
public:
    Normal_Proposer(double rule_scale = HSP_NA_DOUBLE,
                    double rule_from = HSP_NA_DOUBLE,
                    bool is_discrete = false)
        : Proposer("norm", rule_scale, rule_from, {1}, is_discrete)
    {
        user_reset();                                             // initialize
        propose_ = &Normal_Proposer::do_propose_sum_to_previous_; // proposer functional metho
        if (is_rule_start_from_zero())                            // logratio functional method
            logratio_ = &Normal_Proposer::do_logratio_sum_to_zero_;
        else
            logratio_ = &Normal_Proposer::do_logratio_symm_;
    }

    void user_reset() override final
    {
        pDistr_->init_attr("mean", HSP_ATTR_ARGUMENT, 0.);
        pDistr_->init_attr("sd", HSP_ATTR_ARGUMENT, user_rule_scale_);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Lognorm_Proposer : class *** /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Proposer
{
public:
    Lognorm_Proposer(double rule_scale = 1.,
                     double rule_from = HSP_NA_DOUBLE,
                     bool is_discrete = false)
        : Proposer("lnorm", rule_scale, rule_from, {1}, is_discrete)
    {
        user_reset();                                             // initialize
        propose_ = &Lognorm_Proposer::do_propose_times_previous_; // proposer functional method
        if (is_rule_start_from_zero())                            // logratio functional method
            logratio_ = &Lognorm_Proposer::do_logratio_times_zero_;
        else
            logratio_ = &Lognorm_Proposer::do_logratio_lognorm_symm_;
    }

    void user_reset() override final
    {
        pDistr_->init_attr("meanlog", HSP_ATTR_ARGUMENT, 0.);
        pDistr_->init_attr("sdlog", HSP_ATTR_ARGUMENT, user_rule_scale_);
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  proposer switch  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Switchboard to select the custom metropolis proposer rule.
std::shared_ptr<Proposer> proposer_switch(const std::string &rule_name = HSP_NA_STRING,
                                          double rule_scale = HSP_NA_DOUBLE,
                                          double rule_from = HSP_NA_DOUBLE,
                                          bool is_discrete = false)
{
    const str_vec options = {"sample", "unif", "norm", "lnorm"};

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
            "Unknown proposer rule \"" + rule_name + "\"." +
            tools::str::find_closest(rule_name, options, HSP_MIN_STR_DIST));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable : setter  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize variable.
void Variable::init(const std::string &name, double value,
                    bool is_const, bool is_discr,
                    const std::string &func_name, bool is_primitive)
{
    this->name(name);
    user_value_ = HSP_NA_DOUBLE;
    evaluate_ = &Variable::do_eval_error_;

    if (tools::is_na(func_name)) // initialization as a numeric variable
    {
        init_scalar(value);
        user_value_ = value;
    }
    else if (func_name == HSP_IDENTITY_FUNC) // initialization as a pointer
    {
        init_pointer(&user_value_); // NOTE: this is just a placeholder
    }
    else // initialization as a function
    {
        evaluate_ = &Variable::do_eval_func_;
        init_function(compose_func_registries(
            func_name, "function", is_primitive,
            {&UserFuncRegistry::inst(),
             &TranscendentalRegistry::inst(),
             &ArithmeticRegistry::inst(),
             &ContinuosDistrRegistry::inst(),
             &DiscreteDistrRegistry::inst(),
             &InfecProfileRegistry::inst()}));
    }
    is_constant_ = is_const;
    is_discrete_ = is_discr;

    // UPDATE: what if constant = T, discrete = T and value != integer?
}

/// Initialize function reshaping.
/// Drop function if it's constant when reshaping.
void Variable::init_func_reshaping()
{
    if (is_function())
    {
        if (pFunc_->is_constant())
        {
            double x = pFunc_->evaluate();
            evaluate_ = &Variable::do_eval_error_;
            init_scalar(x);
            user_value_ = x;
            is_constant_ = true;
        }
        else
        {
            if (is_discrete_)
                pFunc_->init_discretization();
            pFunc_->init_truncation();

            // FIXME: check that all function arguments have not missing dependencies
            // NOTE: if attributes are NA it could still be fair use as a callable
        }
    }
}

/// Assign distribution with constant arguments.
void Variable::init_distr(const std::string &distr_name)
{
    pDistr_ = compose_distr_registries(
        distr_name,
        "distribution",
        {&ContinuosDistrRegistry::inst(),
         &DiscreteDistrRegistry::inst()});
}

/// Initialize distribution reshaping.
void Variable::init_distr_reshaping()
{
    if (!pDistr_->is_empty())
    {
        if (is_discrete_)
            pDistr_->init_discretization();
        pDistr_->init_truncation();
        // check that non-constant variable has non-zero dimensional range
        if ((!is_constant_) && (pDistr_->get_min() == pDistr_->get_max()))
            throw std::invalid_argument(
                "Degenerate range " + pDistr_->str_range() + " defined for non constant " + err_name() + ".");
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Variable : Setter (Proposer)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign proposer.
void Variable::init_proposer(const std::string &rule_name, double rule_scale, double rule_from)
{
    if (tools::is_na(rule_name)) // if no rule is defined -> constant
    {
        is_constant_ = true;
        propose_ = &Variable::do_propose_error_;
        pProp_ = proposer_switch();
    }
    else // variable to be sampled
    {
        is_constant_ = false;

        if (rule_name == "sample") // sample from model/prior
        {
            pProp_ = std::make_shared<Sample_Proposer>(pDistr_);
            propose_ = &Variable::do_propose_sample_;
            move_logratio_ = &Variable::do_move_logratio_sample_;
            is_valid_move_ = &Variable::do_is_valid_move_sample_;
        }
        else // other proposer rules
        {
            // value from which to start proposal
            if ((!tools::is_na(rule_from)) && (!pDistr_->is_valid(rule_from)))
                tools::warning(
                    "Proposer rule from value " + std::to_string(rule_from) + " out of range" + " {min:" + str_min() + ", max:" + str_max() + "}.");

            // assign the pointer to the proposer
            pProp_ = proposer_switch(rule_name, rule_scale, rule_from, is_discrete());

            // raise a warning if "lnorm" rule when varaible can be negative
            // UPDATE: generalize in case of new positive-defined proposer rules
            if ((pDistr_->get_min() < 0) && (rule_name == "lnorm"))
                tools::warning(
                    "Lognormal proposer rule for non-positive " + err_name() + ".");

            // sampler
            if (is_from_zero())
                propose_ = &Variable::do_propose_from_zero_;
            else
                propose_ = &Variable::do_propose_from_previous_;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable : getter  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Explore the chain of dependencies to check if variable depends already on other one.
/// UPDATE: should also explore the case of functions.
bool Variable::is_chain_dependency(const Variable &other) const
{
    if (name() == other.name())
        return true;
    else if (!is_function())
    {
        if (!other.is_function())
            return get_ptr() == other.get_ptr();
    }
    return false;
}
