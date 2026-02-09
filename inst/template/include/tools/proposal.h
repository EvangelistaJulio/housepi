/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "tools/nan.h"
#include "tools/random.h"
#include "tools/string.h"

#include "function.h"
#include "distribution.h"

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
