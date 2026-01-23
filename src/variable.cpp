/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include "tools/messages.h"
#include "tools/parser.h"
#include "tools/string.h"
#include "tools/random.h"

#include "distribution.h"
#include "proposal.h"
#include "variable.h"

namespace
{
    using std::string;
    using std::vector;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Proposal : constructor  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Empty proposer.
Proposal::Proposal() = default;

/// Proposal constructor from Distribution name.
Proposal::Proposal(const string &distr_name,
                   double proposal_scale,
                   const vector<unsigned short> &args_to_adapt,
                   bool is_discrete)
    : pDistr_(compose_distr_registries(distr_name, "proposal", {"continuos"}))
{
    // set to discrete
    if (is_discrete)
        pDistr_->init_discretization();
    // set range and truncation
    pDistr_->init_truncation();

    // assign user proposal scale
    if (!hsp::na::is(proposal_scale))
        user_scale_ = std::fabs(proposal_scale);

    // load which are the parameter to be adapted
    args_to_adapt_ = args_to_adapt;
}

/// Proposal deconstructor.
Proposal::~Proposal() = default;

/// Return proposer scale.
double Proposal::scale() const { return is_adaptable() ? std::fabs((pDistr_->arg)(args_to_adapt_[0])) : hsp::na::DOUBLE; }

/// Increase spread if acceptation rate is high and viceversa.
void Proposal::adapt(double ratio)
{
    // target accept rate set to 1/4.
    const double factor = 1. + (1. - ratio) * (get_accept_rate() - 0.25);
    // loop over all ids of arguments to adapt
    for (auto i : args_to_adapt_)
        (pDistr_->set_arg)(i, factor * (pDistr_->arg)(i));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  *** Proposal : implementations ***  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  *** Sample_Proposal : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class Sample_Proposal : public Proposal
    {
    public:
        Sample_Proposal(std::shared_ptr<Distribution> distr) { pDistr_ = distr; }

        double propose(hsp::rnd_gen &prng, double) const override final { return pDistr_->rand(prng); }
        double logratio(double xo, double xn) const override final { return pDistr_->logratio(xo, xn); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** Uniform_Proposal : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class Uniform_Proposal : public Proposal
    {
    public:
        Uniform_Proposal(double proposal_scale = hsp::na::DOUBLE, bool is_discrete = false)
            : Proposal("uniform", proposal_scale, {0, 1}, is_discrete) { user_reset(); }

        void user_reset() override final
        {
            pDistr_->init_attr(hsp::params::MIN, Attribute(-user_scale_));
            pDistr_->init_attr(hsp::params::MAX, Attribute(user_scale_));
        }

        double propose(hsp::rnd_gen &prng, double xo) const override final { return xo + pDistr_->rand(prng); }
        double logratio(double, double) const override final { return 0.; }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** Normal_Proposal : class ***  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class Normal_Proposal : public Proposal
    {
    public:
        Normal_Proposal(double proposal_scale = hsp::na::DOUBLE, bool is_discrete = false)
            : Proposal("normal", proposal_scale, {1}, is_discrete) { user_reset(); }

        void user_reset() override final
        {
            pDistr_->init_attr("mean", Attribute(0.));
            pDistr_->init_attr("sd", Attribute(user_scale_));
        }

        double propose(hsp::rnd_gen &prng, double xo) const override final { return xo + pDistr_->rand(prng); }
        double logratio(double, double) const override final { return 0.; }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** Lognorm_Proposal : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class Lognorm_Proposal : public Proposal
    {
    public:
        Lognorm_Proposal(double proposal_scale = 1., bool is_discrete = false)
            : Proposal("lognormal", proposal_scale, {1}, is_discrete) { user_reset(); }

        void user_reset() override final
        {
            pDistr_->init_attr("meanlog", Attribute(0.));
            pDistr_->init_attr("sdlog", Attribute(user_scale_));
        }

        void associate_location(const Attribute &other) override final
        {
            pAttr_loc_ = other.get_ptr();

            if (pAttr_loc_ == nullptr)
            {
                propose_ = &Lognorm_Proposal::do_propose_;
                logratio_ = &Lognorm_Proposal::do_logratio_;
            }
            else
            {
                propose_ = &Lognorm_Proposal::do_propose_loc_;
                logratio_ = &Lognorm_Proposal::do_logratio_loc_;
            }
        };

        double propose(hsp::rnd_gen &prng, double xo) const override final { return (this->*propose_)(prng, xo); }
        double logratio(double xo, double xn) const override final { return (this->*logratio_)(xo, xn); }

    private:
        const Attribute *pAttr_loc_ = nullptr; /// pointer to attribute location

        double (Lognorm_Proposal::*propose_)(hsp::rnd_gen &, double) const = &Lognorm_Proposal::do_propose_;
        double do_propose_(hsp::rnd_gen &prng, double xo) const { return xo * pDistr_->rand(prng); }
        double do_propose_loc_(hsp::rnd_gen &prng, double xo) const
        {
            const double loc = pAttr_loc_->get();
            return loc + (xo - loc) * pDistr_->rand(prng);
        }

        double (Lognorm_Proposal::*logratio_)(double, double) const = &Lognorm_Proposal::do_logratio_;
        double do_logratio_(double xo, double xn) const { return std::log(xo / xn); }
        double do_logratio_loc_(double xo, double xn) const
        {
            const double loc = pAttr_loc_->get();
            return std::log((xo - loc) / (xn - loc));
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  proposer switch  //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Switchboard to select the custom metropolis proposal.
    std::unique_ptr<Proposal> proposer_switch(const string &proposal = "",
                                              double proposal_scale = hsp::na::DOUBLE,
                                              bool is_discrete = false)
    {
        const vector<string> options = {"sample", "uniform", "normal", "lognormal"};

        if (hsp::na::is(proposal))
            return nullptr;

        if (proposal == "uniform")
            return std::make_unique<Uniform_Proposal>(proposal_scale, is_discrete);

        if (proposal == "normal")
            return std::make_unique<Normal_Proposal>(proposal_scale, is_discrete);

        if (proposal == "lognormal")
            return std::make_unique<Lognorm_Proposal>(proposal_scale, is_discrete);

        throw std::invalid_argument(
            "Unknown proposal '" + proposal + "'." + hsp::str::find_closest(proposal, options));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Variable : setter  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Variable::Variable(const std::string &err_type)
{
    this->name("", err_type);
    this->is_constant_ = true;
    buffer_ = hsp::na::DOUBLE;
    user_value_ = hsp::na::DOUBLE;
    value_ = hsp::na::DOUBLE;
    evaluate_ = &Variable::do_eval_error_;
}
Variable::~Variable() = default;

/// Initialize variable.
void Variable::init(double value, bool is_const, std::string_view expr, bool is_primitive)
{
    this->is_constant_ = is_const;
    hsp::parse::Definition parser(expr);

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
void Variable::user_reset(hsp::rnd_gen &prng)
{
    if (!hsp::na::is(user_value_))
        value_ = user_value_;
    else if (is_model())
        init_rand(prng);

    if (pProp_ != nullptr)
        pProp_->user_reset();
}

/// Assign an attribute to function.
void Variable::init_func_attr(std::string_view key, const Attribute &attr) { pFunc_->init_attr(key, attr); }

/// Initialize function reshaping.
/// Drop function if it's constant when reshaping.
/// NOTE: discretization of model has priority
///       e. g. user discr is FALSE, but discrete-defined model
void Variable::init_func_reshaping(bool discr)
{
    if (is_model())
        discr = pDistr_->is_discrete();
    if (discr)
        pFunc_->init_discretization();
    pFunc_->init_truncation();

    // assert no missing arguments
    pFunc_->assert_is_callable_explicitly();

    if (pFunc_->is_constant())
        init(pFunc_->evaluate(), true, "", false);

    else if (is_constant())
        throw std::invalid_argument(
            "Impossible to set " + err_name() + " as constant. Check its dependencies.");

    else if (!pFunc_->is_callable_implicitly())
        get_ = &Variable::do_get_error_;
}

/// Assign distribution with constant arguments.
void Variable::init_distr(std::string_view opt)
{
    pDistr_ = compose_distr_registries(opt, "distribution", {"continuos", "discrete"});
}

/// Assign an attribute to distribution.
void Variable::init_distr_attr(std::string_view key, const Attribute &attr)
{
    pDistr_->init_attr(key, attr);
    if (key == hsp::params::LOC && pProp_ != nullptr)
        pProp_->associate_location(pDistr_->get_loc());
}

/// Initialize distribution reshaping.
void Variable::init_distr_reshaping(bool discr)
{
    if (!is_model())
        return;

    if (discr)
        pDistr_->init_discretization();
    pDistr_->init_truncation();
    pDistr_->assert_is_callable_explicitly();

    // check that non-constant variable has non-zero dimensional range
    helper_degenerate_range_();
}

/// Helper with degenerate range.
void Variable::helper_degenerate_range_()
{
    // skip if one of the range extremes is not provided
    if (pDistr_->get_lower().is_empty() || pDistr_->get_upper().is_empty())
        return;

    // skip if range is not degenerate
    // UPDATE: hard to know in advance if not initialized Variables...
    const auto lower = pDistr_->lower();
    const auto upper = pDistr_->upper();
    if (hsp::na::is(lower) || hsp::na::is(upper) || lower < upper)
        return;

    // case empty interval, lower and upper are "open"
    if (lower > upper || !pDistr_->is_valid(lower) || !pDistr_->is_valid(upper))
        throw std::invalid_argument(
            "Definition produced a degenerate range " + pDistr_->str_range() +
            " for " + err_name() + ".");

    // Warn user if it was not chosen as constant
    if (!is_constant_)
        hsp::msg::warning(
            "Forced at least one instance of " + err_name() +
            " to be constant as its range was degenerate.");

    init_numeric(lower);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Variable : Setter (Proposal)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Random initialization of variable.
void Variable::init_rand(hsp::rnd_gen &prng) { set(pDistr_->rand(prng)); }

/// Assign proposer to variable.
/// NOTE: it must be assigned after model initialization.
void Variable::init_proposer(const string &proposal, double proposal_scale)
{
    is_constant_ = true; // if no proposal is defined, assume constant

    if (!hsp::na::is(proposal))
    {
        is_constant_ = false; // else variable to be sampled

        if (proposal == "sample")
            helper_init_sample_proposer_();
        else
            helper_init_proposer_monte_carlo_(proposal, proposal_scale);
    }
}

double Variable::do_eval_func_(double x) const { return pFunc_->evaluate(x); }

/// Helper with proposer initialization of direct sample from prior/model.
void Variable::helper_init_sample_proposer_()
{
    // assign the pointer to the proposer
    pProp_ = std::make_unique<Sample_Proposal>(pDistr_);

    // setup sample proposal functional methods
    get_move_logratio_ = &Variable::do_move_logratio_sample_;
    is_valid_move_ = &Variable::do_is_valid_move_always_;
}

/// Helper with Monte Carlo proposer initialization.
void Variable::helper_init_proposer_monte_carlo_(const string &proposal, double proposal_scale)
{
    // raise a warning if "lnorm" proposal when variable can be negative
    if (pDistr_ != nullptr && pDistr_->lower() < 0 && proposal == "lnorm")
        hsp::msg::warning(
            "Provided a lognormal proposal for non-positive " + err_name() + ".");

    // select proposer
    pProp_ = proposer_switch(proposal, proposal_scale, is_discrete());

    // setup parameteric proposal functional methods
    get_move_logratio_ = &Variable::do_get_move_logratio_default_;
    if (pDistr_ == nullptr)
        is_valid_move_ = &Variable::do_is_valid_move_always_;
    else
        is_valid_move_ = &Variable::do_is_valid_move_default_;
}

/// Wheter the variable is constant.
bool Variable::is_constant() const noexcept
{
    if (is_numeric())
        return is_constant_;
    if (is_pointer())
        return false;
    if (is_function())
        return pFunc_->is_constant();
    return is_constant_;
}
/// Wheter the variable is discrete.
bool Variable::is_discrete() const noexcept
{
    if (is_model())
        return pDistr_->is_discrete();
    return Attribute::is_discrete();
}
/// Wheter the variable has a prior or a model associated.
bool Variable::is_model() const noexcept { return pDistr_ != nullptr && !pDistr_->is_empty(); }
/// Wheter the variable has an adaptable parametric proposal rule.
bool Variable::is_adaptable() const noexcept { return pProp_ != nullptr && pProp_->is_adaptable(); }

/// Set new proposal internally and store a buffer of old value.
void Variable::propose(hsp::rnd_gen &prng) { set(pProp_->propose(prng, value_)); }

double Variable::get_distr_logpdf() const { return pDistr_->logpdf(value_); }
double Variable::get_distr_logratio() const { return pDistr_->logratio(buffer_, value_); }
double Variable::get_propose_logratio() const { return pProp_->logratio(buffer_, value_); }

void Variable::update_accept_rate(bool is_accepted) { pProp_->update_accept_rate(is_accepted); }
void Variable::reset_accept_rate() { pProp_->reset_accept_rate(); }
void Variable::adapt_proposal(double ratio) { pProp_->adapt(ratio); }

/// Get a description about the functional status.
string Variable::str_summary() const
{
    auto output = Attribute::str_summary();
    if (is_model())
        output += " | ~" + pDistr_->str_summary(false);
    return output;
}

bool Variable::do_is_valid_move_default_() const { return pDistr_->is_valid(value_); }
