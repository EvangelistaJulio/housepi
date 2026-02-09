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
#include "tools/nan.h"
#include "tools/parser.h"
#include "tools/proposal.h"
#include "tools/string.h"
#include "tools/random.h"

#include "variable.h"
#include "function.h"
#include "distribution.h"

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
