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

    // load which are the parameter to be adapted
    n_args_to_adapt_ = args_to_adapt.size();
    args_to_adapt_list_ = new unsigned short[n_args_to_adapt_];

    for (unsigned short c = 0; c < n_args_to_adapt_; ++c)
        args_to_adapt_list_[c] = args_to_adapt[c];
};

/// Proposer constructor from Distribution name.
Proposer::Proposer(const std::string &name,
                   double rule_scale,
                   const std::vector<unsigned short> &args_to_adapt,
                   bool is_discrete)
    : Proposer(ContinuosDistrRegistry::inst().get_distr(name),
               rule_scale, args_to_adapt, is_discrete) {}

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

// NOTE: this is just a place holder within Variable.
class Sample_Proposer : public Proposer
{
public:
    Sample_Proposer() { name_ = "sample"; }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Uniform_Proposer : class *** /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Proposer
{
public:
    Uniform_Proposer(double rule_scale = HSP_NA_DOUBLE, bool is_discrete = false)
        : Proposer("unif", rule_scale, {0, 1}, is_discrete) { user_reset(); }

    void user_reset() override final
    {
        pDistr_->init_attr(HSP_PARAM_MIN, -user_rule_scale_);
        pDistr_->init_attr(HSP_PARAM_MAX, user_rule_scale_);
    }

    double propose(rnd_gen &prng, double xo) override final { return xo + pDistr_->rand(prng); }
    double logratio(double, double) const override final { return 0.; }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Normal_Proposer : class *** //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Proposer
{
public:
    Normal_Proposer(double rule_scale = HSP_NA_DOUBLE, bool is_discrete = false)
        : Proposer("norm", rule_scale, {1}, is_discrete) { user_reset(); }

    void user_reset() override final
    {
        pDistr_->init_attr("mean", 0.);
        pDistr_->init_attr("sd", user_rule_scale_);
    }

    double propose(rnd_gen &prng, double xo) override final { return xo + pDistr_->rand(prng); }
    double logratio(double, double) const override final { return 0.; }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Lognorm_Proposer : class *** /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Proposer
{
public:
    Lognorm_Proposer(double rule_scale = 1., bool is_discrete = false)
        : Proposer("lnorm", rule_scale, {1}, is_discrete) { user_reset(); }

    void user_reset() override final
    {
        pDistr_->init_attr("meanlog", 0.);
        pDistr_->init_attr("sdlog", user_rule_scale_);
    }

    double propose(rnd_gen &prng, double xo) override final { return xo * pDistr_->rand(prng); }
    double logratio(double xo, double xn) const override final { return std::log(xo / xn); }
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  proposer switch  //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Switchboard to select the custom metropolis proposer rule.
std::shared_ptr<Proposer> proposer_switch(const std::string &rule_name = HSP_NA_STRING,
                                          double rule_scale = HSP_NA_DOUBLE,
                                          bool is_discrete = false)
{
    const str_vec options = {"sample", "unif", "norm", "lnorm"};

    if (tools::is_na(rule_name))
        return std::make_shared<Proposer>(std::make_shared<Distribution>());

    else if (rule_name == "unif")
        return std::make_shared<Uniform_Proposer>(rule_scale, is_discrete);

    else if (rule_name == "norm")
        return std::make_shared<Normal_Proposer>(rule_scale, is_discrete);

    else if (rule_name == "lnorm")
        return std::make_shared<Lognorm_Proposer>(rule_scale, is_discrete);

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
                    const DefinitionParser &parser, bool is_primitive)
{
    this->name_ = name;
    this->is_constant_ = is_const;
    this->is_discrete_ = is_discr;
    // UPDATE: what if constant = T, discrete = T and value != integer?

    user_value_ = HSP_NA_DOUBLE;
    evaluate_ = &Variable::do_eval_error_;

    if (parser.is_na()) // initialization as a numeric variable
    {
        init_scalar(value);
        user_value_ = value;
    }
    else if (parser.is_terminal()) // initialization as a pointer
    {
        init_pointer(&user_value_); // NOTE: this is just a placeholder
    }
    else // initialization as a function
    {
        init_function(parser.root.func, is_primitive);
        evaluate_ = &Variable::do_eval_func_;
    }
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
void Variable::init_distr(const DefinitionParser &expr)
{
    pDistr_ = compose_distr_registries(
        expr.root.func, 
        "distribution",
        {&ContinuosDistrRegistry::inst(),
         &DiscreteDistrRegistry::inst()});
}
/// Assign an attribute to distribution.
void Variable::init_distr_attr(const std::string &key, const Attribute &attr)
{
    if (key == HSP_PARAM_LOC)
        init_location_arg(attr);
    else
        pDistr_->init_attr(key, attr);
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
                "Degenerate range " + pDistr_->str_range() +
                " defined for non constant " + err_name() + ".");
    }
}

/// Random initialization of variable.
void Variable::init_rand(rnd_gen &prng)
{
    if (loc_.is_empty())
        set(do_propose_sample_(prng));
    else
        set(do_propose_sample_with_loc_(prng));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Variable : Setter (Proposer)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize location parameter.
void Variable::init_location_arg(const Attribute &x)
{
    this->loc_ = x;

    // assign location parameter
    if (loc_.is_empty())
    {
        distr_logpdf_ = &Variable::do_distr_logpdf_default_;
        distr_logratio_ = &Variable::do_distr_logratio_default_;
    }
    else
    {
        distr_logpdf_ = &Variable::do_distr_logpdf_with_loc_;
        distr_logratio_ = &Variable::do_distr_logratio_with_loc_;
    }
}

/// Assign proposer.
void Variable::init_proposer(const std::string &rule_name, double rule_scale)
{
    if (tools::is_na(rule_name))
    {
        // if no rule is defined assume constant
        is_constant_ = true;
        propose_ = &Variable::do_propose_error_;
        pProp_ = proposer_switch();
    }
    else
    {
        // variable to be sampled
        is_constant_ = false;

        if (rule_name == "sample")
            helper_init_sample_proposer_();
        else // other proposer rules
        {
            // raise a warning if "lnorm" rule when varaible can be negative
            // UPDATE: generalize in case of new positive-defined proposer rules
            if ((pDistr_->get_min() < 0) && (rule_name == "lnorm"))
                tools::warning(
                    "Lognormal proposer rule for non-positive " + err_name() + ".");

            // assign the pointer to the proposer
            pProp_ = proposer_switch(rule_name, rule_scale, is_discrete());

            if (loc_.is_empty())
            {
                propose_ = &Variable::do_propose_default_;
                is_valid_move_ = &Variable::do_is_valid_move_default_;
                get_propose_logratio_ = &Variable::do_get_propose_logratio_default_;
                get_move_logratio_ = &Variable::do_get_move_logratio_default_;
            }
            else
            {
                propose_ = &Variable::do_propose_with_loc_;
                is_valid_move_ = &Variable::do_is_valid_move_with_loc_;
                get_propose_logratio_ = &Variable::do_get_propose_logratio_with_loc_;
                get_move_logratio_ = &Variable::do_get_move_logratio_with_loc_;
            }
        }
    }
}

/// Helper with proposer initialization of direct sample from prior/model.
void Variable::helper_init_sample_proposer_()
{
    pProp_ = std::make_shared<Sample_Proposer>();

    if (loc_.is_empty())
    {
        propose_ = &Variable::do_propose_sample_;
        get_propose_logratio_ = &Variable::do_get_propose_logratio_sample_;
    }
    else
    {
        propose_ = &Variable::do_propose_sample_with_loc_;
        get_propose_logratio_ = &Variable::do_get_propose_logratio_sample_with_loc_;
    }

    is_valid_move_ = &Variable::do_is_valid_move_sample_;
    get_move_logratio_ = &Variable::do_move_logratio_sample_;
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
    else if (!(is_function() || other.is_function()))
        return get_ptr() == other.get_ptr();
    else
        return false;
}
