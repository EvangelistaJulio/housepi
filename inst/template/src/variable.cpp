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
#include "tools/error.h"
#include "tools/params.h"
#include "tools/string.h"
#include "distribution_registry.h"

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
    if (!hsp::na::is(rule_scale))
        user_rule_scale_ = std::fabs(rule_scale);

    // load which are the parameter to be adapted
    args_to_adapt_ = args_to_adapt;
};

/// Proposer constructor from Distribution name.
Proposer::Proposer(const std::string &distr_name,
                   double rule_scale,
                   const std::vector<unsigned short> &args_to_adapt,
                   bool is_discrete)
    : Proposer(ContinuosDistrRegistry::inst().get_distr(distr_name),
               rule_scale, args_to_adapt, is_discrete) {}

/// Increase spread if acceptation rate is high and viceversa.
void Proposer::adapt(double ratio)
{
    // target accept rate set to 1/4.
    const double factor = 1. + (1. - ratio) * (get_accept_rate() - 0.25);
    // loop over all ids of arguments to adapt
    for (unsigned short i : args_to_adapt_)
        (pDistr_->set_arg)(i, factor * (pDistr_->arg)(i));
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
    Uniform_Proposer(double rule_scale = hsp::na::DOUBLE, bool is_discrete = false)
        : Proposer("uniform", rule_scale, {0, 1}, is_discrete) { user_reset(); }

    void user_reset() override final
    {
        pDistr_->init_attr(hsp::params::MIN, Attribute(-user_rule_scale_));
        pDistr_->init_attr(hsp::params::MAX, Attribute(user_rule_scale_));
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
    Normal_Proposer(double rule_scale = hsp::na::DOUBLE, bool is_discrete = false)
        : Proposer("normal", rule_scale, {1}, is_discrete) { user_reset(); }

    void user_reset() override final
    {
        pDistr_->init_attr("mean", Attribute(0.));
        pDistr_->init_attr("sd", Attribute(user_rule_scale_));
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
        : Proposer("lognormal", rule_scale, {1}, is_discrete) { user_reset(); }

    void user_reset() override final
    {
        pDistr_->init_attr("meanlog", Attribute(0.));
        pDistr_->init_attr("sdlog", Attribute(user_rule_scale_));
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
std::shared_ptr<Proposer> proposer_switch(const std::string &rule_name = hsp::na::STRING,
                                          double rule_scale = hsp::na::DOUBLE,
                                          bool is_discrete = false)
{
    const str_vec options = {"sample", "uniform", "normal", "lognormal"};

    if (hsp::na::is(rule_name))
        return std::make_shared<Proposer>(std::make_shared<Distribution>());

    else if (rule_name == "uniform")
        return std::make_shared<Uniform_Proposer>(rule_scale, is_discrete);

    else if (rule_name == "normal")
        return std::make_shared<Normal_Proposer>(rule_scale, is_discrete);

    else if (rule_name == "lognormal")
        return std::make_shared<Lognorm_Proposer>(rule_scale, is_discrete);

    throw std::invalid_argument(
        "Unknown proposer rule \"" + rule_name + "\"." +
        hsp::str::find_closest(rule_name, options));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable : setter  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize variable.
void Variable::init(double value, bool is_const,
                    const DefinitionParser &parser, bool is_primitive)
{
    this->is_constant_ = is_const;

    user_value_ = hsp::na::DOUBLE;
    evaluate_ = &Variable::do_eval_error_;

    if (is_const || parser.is_na()) // initialization as a numeric variable
    {
        init_numeric(value);
        user_value_ = value;
        // UPDATE: what if constant = T, discrete = T and value != integer?
    }
    else if (parser.is_terminal()) // initialization as a pointer
    {
        init_pointer(this); // NOTE: this is just a placeholder
    }
    else // initialization as a function
    {
        init_function(parser.root.func, is_primitive);
        evaluate_ = &Variable::do_eval_func_;
    }
}

/// Restore user defined variable initialization.
///
/// NOTE: by construction convention it is possible to randomly initialize
///       if user value is NA but a distribution is given.
///       Otherwise it is truly a constant NA value.
void Variable::user_reset(rnd_gen &prng)
{
    if (!hsp::na::is(user_value_))
        value_ = user_value_;
    else if (!pDistr_->is_empty())
        init_rand(prng);

    pProp_->user_reset();
}

/// Initialize function reshaping.
/// Drop function if it's constant when reshaping.
void Variable::init_func_reshaping(bool discr)
{
    // NOTE: discretization of model has priority
    // e. g. user discr is FALSE, but discrete-defined model
    if (!pDistr_->is_empty())
        discr = pDistr_->is_discrete();
    if (discr)
        pFunc_->init_discretization();
    pFunc_->init_truncation();

    if (pFunc_->is_constant())
        init(pFunc_->evaluate(), true);

    else if (is_constant())
        throw std::invalid_argument(
            "Impossible to set " + err_name() + " as constant. Check its dependencies.");

    else if (!pFunc_->is_implicit_callable())
        get_ = &Variable::do_get_error_;
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
    if (key == hsp::params::LOCATION)
        init_location_arg(attr);
    else
        pDistr_->init_attr(key, attr);
}

/// Initialize distribution reshaping.
void Variable::init_distr_reshaping(bool discr)
{
    if (pDistr_->is_empty())
        return;

    if (discr)
        pDistr_->init_discretization();
    pDistr_->init_truncation();

    // check that non-constant variable has non-zero dimensional range
    if (!is_constant_ && pDistr_->get_min() == pDistr_->get_max())
        throw std::invalid_argument(
            "Degenerate range " + pDistr_->str_range() +
            " defined for non constant " + err_name() + ".");
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

/// Assign proposer to variable.
/// NOTE: it must be assigned after model initialization.
void Variable::init_proposer(const std::string &rule_name, double rule_scale)
{
    if (hsp::na::is(rule_name))
    {
        // if no rule is defined assume constant
        is_constant_ = true;
        propose_ = &Variable::do_propose_error_;
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
            if (pDistr_->get_min() < 0 && rule_name == "lnorm")
                hsp::msg::warning(
                    "Provided a lognormal proposer rule for non-positive " + err_name() + ".");

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

/// Wheter the variable is constant.
bool Variable::is_constant() const
{
    if (is_numeric())
        return is_constant_;
    else if (is_pointer())
        return false;
    else if (is_function())
        return pFunc_->is_constant();
    else
        return is_constant_;
}
/// Wheter the variable is discrete.
bool Variable::is_discrete() const
{
    if (pDistr_->is_empty())
        return Attribute::is_discrete();
    else
        return pDistr_->is_discrete();
}

/// Get a description about the functional status.
std::string Variable::str_summary() const
{
    std::string output = Attribute::str_summary();
    if (is_model())
        output += " | ~" + pDistr_->str_summary(false);
    return output;
}
