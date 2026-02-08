/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <string>
#include <vector>

#include "tools/messages.h"
#include "tools/parser.h"
#include "tools/string.h"
#include "tools/random.h"

#include "attribute.h"
#include "distribution.h"
#include "variable.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Proposal : class ***  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Base class for Metropolis proposals.
class Proposal
{
public:
    Proposal() = default;
    virtual ~Proposal() = default;

    virtual void user_reset() {}
    virtual void associate_location(const Attribute &) {}

    virtual double propose(hsp::rnd_gen &, double) const = 0;
    virtual double logratio(double, double) const = 0;

    virtual void update_accept_rate(bool) noexcept {};
    virtual void reset_accept_rate() noexcept {};
    virtual void adapt(double) {};

    virtual bool is_adaptable() const noexcept = 0;
    virtual std::string_view name() const noexcept = 0;
};

namespace
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  *** Sample_Proposal : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class Sample_Proposal final : public Proposal
    {
    public:
        Sample_Proposal(const Distribution *distr) : pDistr_(distr) {}

        double propose(hsp::rnd_gen &prng, double) const override final { return pDistr_->rand(prng); }
        double logratio(double xo, double xn) const override final { return pDistr_->logratio(xo, xn); }

        bool is_adaptable() const noexcept override final { return false; }
        std::string_view name() const noexcept override final { return "sample"; }

    private:
        const Distribution *pDistr_;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** AdaptableProposal : class ***  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Adaptable version of Metropolis proposals.
    class AdaptableProposal : public Proposal
    {
    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        AdaptableProposal(const std::string &distr_name, double scale, bool is_discrete)
            : pDistr_(compose_distr_registries(distr_name, "proposal", {"continuos"}))
        {
            if (is_discrete)
                pDistr_->init_discretization();
            pDistr_->init_truncation();

            if (!hsp::na::is(scale))
                user_scale_ = std::fabs(scale);
        }
        virtual ~AdaptableProposal() = default;

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        void update_accept_rate(bool is_accepted) noexcept override final
        {
            accept_rate_ = static_cast<double>(accept_rate_ * n_proposed_ + is_accepted) / (1. + n_proposed_);
            ++n_proposed_;
        }
        void reset_accept_rate() noexcept override final
        {
            accept_rate_ = 0.;
            n_proposed_ = 0;
        }

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Wheter proposer can be adapted.
        bool is_adaptable() const noexcept override final { return true; }

        ////////////////////
        //////  SHOW  //////
        ////////////////////

        std::string_view name() const noexcept override final { return pDistr_->name(); }

    protected:
        std::unique_ptr<Distribution> pDistr_; /// distribution associated to proposal
        double user_scale_ = 1.;               /// scale of the Metropolis proposal
        double accept_rate_ = 0.;              /// acceptance rate
        unsigned long n_proposed_ = 0;         /// n. of proposed moves

        /// Factor used for adaptation, with target accept rate set to 25%.
        double factor_(double ratio) { return 1. + (1. - ratio) * (accept_rate_ - 0.25); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** Uniform_Proposal : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class Uniform_Proposal final : public AdaptableProposal
    {
    public:
        Uniform_Proposal(double scale, bool is_discrete)
            : AdaptableProposal("uniform", scale, is_discrete) { user_reset(); }

        void user_reset() override final
        {
            pDistr_->init_attr("min", -user_scale_);
            pDistr_->init_attr("max", user_scale_);
        }

        double propose(hsp::rnd_gen &prng, double xo) const override final { return xo + pDistr_->rand(prng); }
        double logratio(double, double) const override final { return 0.; }

        void adapt(double ratio) override final
        {
            const double m = factor_(ratio);
            pDistr_->set_lower(m * pDistr_->lower().value());
            pDistr_->set_upper(m * pDistr_->upper().value());
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** Normal_Proposal : class ***  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class Normal_Proposal final : public AdaptableProposal
    {
    public:
        Normal_Proposal(double scale, bool is_discrete)
            : AdaptableProposal("normal", scale, is_discrete) { user_reset(); }

        void user_reset() override final
        {
            pDistr_->init_attr("mean", 0.);
            pDistr_->init_attr("sd", user_scale_);
        }

        double propose(hsp::rnd_gen &prng, double xo) const override final { return xo + pDistr_->rand(prng); }
        double logratio(double, double) const override final { return 0.; }

        /// Adapt standard deviation.
        void adapt(double ratio) override final { pDistr_->set_arg(1, factor_(ratio) * pDistr_->arg(1)); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** Lognorm_Proposal : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class Lognorm_Proposal final : public AdaptableProposal
    {
    public:
        Lognorm_Proposal(double scale, bool is_discrete)
            : AdaptableProposal("lognormal", scale, is_discrete) { user_reset(); }

        void user_reset() override final
        {
            pDistr_->init_attr("meanlog", 0.);
            pDistr_->init_attr("sdlog", user_scale_);
        }

        void associate_location(const Attribute &other) override final
        {
            pAttr_loc_ = other.pointer();

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

        /// Adapt log of standard deviation.
        void adapt(double ratio) override final { pDistr_->set_arg(1, factor_(ratio) * pDistr_->arg(1)); }

    private:
        const Attribute *pAttr_loc_ = nullptr; /// pointer to attribute location

        double (Lognorm_Proposal::*propose_)(hsp::rnd_gen &, double) const = &Lognorm_Proposal::do_propose_;
        double do_propose_(hsp::rnd_gen &prng, double xo) const { return xo * pDistr_->rand(prng); }
        double do_propose_loc_(hsp::rnd_gen &prng, double xo) const
        {
            const double loc = pAttr_loc_->value();
            return loc + (xo - loc) * pDistr_->rand(prng);
        }

        double (Lognorm_Proposal::*logratio_)(double, double) const = &Lognorm_Proposal::do_logratio_;
        double do_logratio_(double xo, double xn) const { return std::log(xo / xn); }
        double do_logratio_loc_(double xo, double xn) const
        {
            const double loc = pAttr_loc_->value();
            return std::log((xo - loc) / (xn - loc));
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////  proposer switch  //////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Switchboard to select the custom metropolis proposal.
    std::shared_ptr<Proposal> proposal_switch(const std::string &name, double scale, bool is_discrete)
    {
        const std::vector<std::string> options = {"sample", "uniform", "normal", "lognormal"};

        if (hsp::na::is(name))
            return nullptr;

        if (name == "uniform")
            return std::make_shared<Uniform_Proposal>(scale, is_discrete);

        if (name == "normal")
            return std::make_shared<Normal_Proposal>(scale, is_discrete);

        if (name == "lognormal")
            return std::make_shared<Lognorm_Proposal>(scale, is_discrete);

        throw std::invalid_argument(
            "Unknown proposal '" + name + "'." + hsp::str::find_closest(name, options));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Variable : class ***  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Variable::~Variable() = default;

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
    if (key == "loc" && pProp_ != nullptr)
        pProp_->associate_location(pDistr_->loc());
}

/// Initialize distribution reshaping.
void Variable::init_distr_reshaping(bool is_discr)
{
    if (!is_model())
        return;

    if (is_discr)
        pDistr_->init_discretization();
    pDistr_->init_truncation();

    pDistr_->assert_is_callable_explicitly();

    // check that non-constant variable has non-zero dimensional range
    helper_assert_non_degenerate_range_();
}

/// Helper with degenerate range.
void Variable::helper_assert_non_degenerate_range_()
{
    // skip if one of the range extremes is not provided
    if (!pDistr_->lower().is_callable() || !pDistr_->upper().is_callable())
        return;

    // skip if range is not degenerate
    const double sx = pDistr_->lower().value();
    const double dx = pDistr_->upper().value();
    if (sx < dx)
        return;

    // case empty interval, lower and upper are "open"
    if (sx > dx || !pDistr_->is_valid(sx) || !pDistr_->is_valid(dx))
        throw std::invalid_argument(
            "Definition produced a inconsistent range " + pDistr_->str_range() +
            " for " + err_name() + ".");

    if (!pDistr_->lower().is_constant() || !pDistr_->upper().is_constant())
        return;

    // warn user if it was not chosen as constant
    if (!is_user_const_)
    {
        hsp::msg::warning(
            "Turned at least one " + err_name() + " to constant as range was degenerate.");
        init_numeric(sx);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Variable : Setter (Proposal)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign proposer to variable.
/// NOTE: it must be assigned after model initialization.
void Variable::init_proposal(const std::string &proposal, double proposal_scale)
{
    if (mode_ != Mode::Empty && mode_ != Mode::Numeric)
        return;

    if (!hsp::na::is(proposal))
    {
        if (mode_ == Mode::Empty)
            init_numeric(hsp::na::DOUBLE);
        is_user_const_ = false;

        if (proposal == "sample")
            helper_init_sample_proposal_();
        else
            helper_init_monte_carlo_proposal_(proposal, proposal_scale);
    }
}

/// Helper with proposer initialization of direct sample from prior/model.
void Variable::helper_init_sample_proposal_()
{
    // assign the pointer to the proposer
    pProp_ = std::make_shared<Sample_Proposal>(pDistr_.get());

    // setup sample proposal functional methods
    ratio_mode_ = Ratio::None;
    valid_mode_ = Validity::Always;
}

/// Helper with Monte Carlo proposer initialization.
void Variable::helper_init_monte_carlo_proposal_(const std::string &proposal, double proposal_scale)
{
    // raise a warning if "lnorm" proposal when variable can be negative
    if (pDistr_ != nullptr && pDistr_->lower().value() < 0 && proposal == "lnorm")
        hsp::msg::warning(
            "Provided a lognormal proposal for non-positive " + err_name() + ".");

    // select proposer
    pProp_ = proposal_switch(proposal, proposal_scale, is_discrete());

    // set modes
    ratio_mode_ = Ratio::Compute;
    valid_mode_ = pDistr_ == nullptr ? Validity::Always : Validity::Check;
}

/// Restore user defined variable initialization.
void Variable::user_reset(hsp::rnd_gen &prng)
{
    if (is_constant())
        return;

    if (is_model())
        set(pDistr_->rand(prng));

    if (pProp_ != nullptr)
        pProp_->user_reset();
}

void Variable::optimize()
{
    if (is_constant())
        init_numeric(value());
    else if (is_user_const_)
        throw std::invalid_argument(
            "Impossible to set " + err_name() + " as constant. Check its dependencies.");

    // reduce function indirection
    if (is_pointer() && pointer()->is_function())
        init_function(pointer()->function());
}

/// Sugar to call variable as a function if applicable, otherwise raise error.
double Variable::operator()(double x) const
{
    switch (mode_)
    {
    case Mode::Empty:
    case Mode::Numeric:
    case Mode::Pointer:
        throw std::runtime_error(
            "Called undefined evaluate operator for non-function " + err_name() + ".");
    case Mode::Function:
        return pFunc_->evaluate(x);
    }
    __builtin_unreachable();
}

bool Variable::is_valid_move() const
{
    switch (valid_mode_)
    {
    case Validity::Always:
        return true;
    case Validity::Check:
        return pDistr_->is_valid(value_);
    }
    __builtin_unreachable();
}

/// Wheter the variable is discrete.
bool Variable::is_discrete() const noexcept { return is_model() ? pDistr_->is_discrete() : Attribute::is_discrete(); }
/// Wheter the variable has a prior or a model associated.
bool Variable::is_model() const noexcept { return pDistr_ != nullptr; }
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
std::string Variable::str_summary() const
{
    return Attribute::str_summary() + (is_model() ? " | ~" + pDistr_->str_summary(false) : "");
}
