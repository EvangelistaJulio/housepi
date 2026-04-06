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
#include "tools/string.h"
#include "tools/random.h"

#include "distribution.h"
#include "function.h"
#include "registry.h"
#include "variable.h"

namespace hsp
{
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** Proposal : class ***  /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Base class for Metropolis proposals.
    class Proposal
    {
    public:
        /// Default constructor.
        Proposal() = default;
        /// Virtual destructor.
        virtual ~Proposal() = default;
        /// Clone this proposal object.
        virtual std::unique_ptr<Proposal> clone() const = 0;

        /// Reset to user-specified scale.
        virtual void user_reset() {}
        /// Associate proposal with a location attribute.
        virtual void associate_location(const Attribute &) {}

        /// Propose a value using random number generator.
        virtual double propose(rnd_gen &, double) const = 0;
        /// Compute log-ratio of proposal densities for Metropolis-Hastings.
        virtual double log_ratio(double, double) const = 0;

        /// Update acceptance rate with result of last proposal.
        virtual void update_accept_rate(bool) noexcept {};
        /// Reset acceptance rate counter.
        virtual void reset_accept_rate() noexcept {};
        /// Adapt proposal parameters based on acceptance ratio.
        virtual void adapt(double) {};

        /// Check if proposal is adaptable.
        virtual bool is_adaptable() const noexcept = 0;
        /// Get name of proposal distribution.
        virtual std::string_view name() const noexcept = 0;

        /// Rebind proposal to a new distribution.
        virtual void rebind_distrib(const Distribution *) {}
    };

#define MAKE_PROPOSAL_CLONE_METHOD(Class_) \
    std::unique_ptr<Proposal> clone() const override { return std::make_unique<Class_>(*this); }
#define MAKE_PROPOSAL_NAME_METHOD(name_) \
    std::string_view name() const noexcept override final { return #name_; }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** NullProposal_ : class ***  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class NullProposal_ final : public Proposal
    {
    public:
        NullProposal_() {}
        MAKE_PROPOSAL_CLONE_METHOD(NullProposal_)
        MAKE_PROPOSAL_NAME_METHOD()

        double propose(rnd_gen &, double xo) const override final { return xo; }
        double log_ratio(double, double) const override final { return 0.; }

        bool is_adaptable() const noexcept override final { return false; }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  *** SampleProposal_ : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class SampleProposal_ final : public Proposal
    {
    public:
        SampleProposal_(const Distribution *distr) : pDistr_(distr) {}
        MAKE_PROPOSAL_CLONE_METHOD(SampleProposal_)
        MAKE_PROPOSAL_NAME_METHOD(sample)

        void rebind_distrib(const Distribution *new_distr) override final { pDistr_ = new_distr; }

        double propose(rnd_gen &prng, double) const override final { return pDistr_->rand(prng); }
        double log_ratio(double xo, double xn) const override final { return pDistr_->log_ratio(xo, xn); }

        bool is_adaptable() const noexcept override final { return false; }

    private:
        const Distribution *pDistr_;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** AdaptableProposal : class ***  ///////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Threshold for proposal scale.
    constexpr double PROPOSAL_XTOL_SCALE = 1.e-8;

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

            if (!na::is(scale))
                user_scale_ = std::fabs(scale);
        }
        virtual ~AdaptableProposal() = default;

        AdaptableProposal(const AdaptableProposal &other)
            : pDistr_(other.pDistr_ ? other.pDistr_->clone_distrib() : nullptr),
              user_scale_(other.user_scale_),
              accept_rate_(other.accept_rate_),
              n_proposed_(other.n_proposed_) {}

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

        /// Adapt (positive) scale with target accept rate set to 25%.
        /// Ref: Andrieu and Thomas (2008)
        double adapt_scale_(double xold, double ratio)
        {
            xold *= 1. + (1. - ratio) * (accept_rate_ - 0.25);
            return (xold < PROPOSAL_XTOL_SCALE) ? PROPOSAL_XTOL_SCALE : xold;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** UniformProposal_ : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class UniformProposal_ final : public AdaptableProposal
    {
    public:
        UniformProposal_(double scale, bool is_discrete)
            : AdaptableProposal("uniform", scale, is_discrete) { user_reset(); }
        MAKE_PROPOSAL_CLONE_METHOD(UniformProposal_)

        void user_reset() override final
        {
            pDistr_->init_attr("min", -user_scale_);
            pDistr_->init_attr("max", user_scale_);
        }

        double propose(rnd_gen &prng, double xo) const override final { return xo + pDistr_->rand(prng); }
        double log_ratio(double, double) const override final { return 0.; }

        void adapt(double ratio) override final
        {
            const double m = adapt_scale_(pDistr_->upper().value(), ratio);
            pDistr_->set_lower(-m);
            pDistr_->set_upper(m);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** NormalProposal_ : class ***  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class NormalProposal_ final : public AdaptableProposal
    {
    public:
        NormalProposal_(double scale, bool is_discrete)
            : AdaptableProposal("normal", scale, is_discrete) { user_reset(); }
        MAKE_PROPOSAL_CLONE_METHOD(NormalProposal_)

        void user_reset() override final
        {
            pDistr_->init_attr("mean", 0.);
            pDistr_->init_attr("sd", user_scale_);
        }

        double propose(rnd_gen &prng, double xo) const override final { return xo + pDistr_->rand(prng); }
        double log_ratio(double, double) const override final { return 0.; }

        /// Adapt standard deviation.
        void adapt(double ratio) override final { pDistr_->set_arg(1, adapt_scale_(pDistr_->arg(1), ratio)); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** LognormProposal_ : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class LognormProposal_ final : public AdaptableProposal
    {
    public:
        LognormProposal_(double scale, bool is_discrete)
            : AdaptableProposal("lognormal", scale, is_discrete) { user_reset(); }
        MAKE_PROPOSAL_CLONE_METHOD(LognormProposal_)

        void user_reset() override final
        {
            pDistr_->init_attr("meanlog", 0.);
            pDistr_->init_attr("sdlog", user_scale_);
        }

        void associate_location(const Attribute &other) override final
        {
            pAttr_loc_ = other.pointer();

            if (pAttr_loc_)
            {
                propose_ = &LognormProposal_::do_propose_loc_;
                logratio_ = &LognormProposal_::do_logratio_loc_;
            }
            else
            {
                propose_ = &LognormProposal_::do_propose_;
                logratio_ = &LognormProposal_::do_logratio_;
            }
        };

        double propose(rnd_gen &prng, double xo) const override final { return (this->*propose_)(prng, xo); }
        double log_ratio(double xo, double xn) const override final { return (this->*logratio_)(xo, xn); }

        /// Adapt log of standard deviation.
        void adapt(double ratio) override final { pDistr_->set_arg(1, adapt_scale_(pDistr_->arg(1), ratio)); }

    private:
        const Attribute *pAttr_loc_ = nullptr; /// pointer to attribute location

        double (LognormProposal_::*propose_)(rnd_gen &, double) const = &LognormProposal_::do_propose_;
        double do_propose_(rnd_gen &prng, double xo) const { return xo * pDistr_->rand(prng); }
        double do_propose_loc_(rnd_gen &prng, double xo) const
        {
            const double loc = pAttr_loc_->value();
            return loc + (xo - loc) * pDistr_->rand(prng);
        }

        double (LognormProposal_::*logratio_)(double, double) const = &LognormProposal_::do_logratio_;
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
    std::unique_ptr<Proposal> proposal_switch(const std::string &name, double scale, bool is_discrete)
    {
        const std::vector<std::string> options = {"sample", "uniform", "normal", "lognormal"};

        if (name == "uniform")
            return std::make_unique<UniformProposal_>(scale, is_discrete);

        if (name == "normal")
            return std::make_unique<NormalProposal_>(scale, is_discrete);

        if (name == "lognormal")
            return std::make_unique<LognormProposal_>(scale, is_discrete);

        throw std::invalid_argument(
            "Unknown proposal '" + name + "'." + str::find_closest(name, options));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** Variable : class ***  /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    Variable::Variable() = default;
    Variable::~Variable() = default;

    Variable::Variable(const Variable &other)
        : Attribute(other),
          valid_mode_(other.valid_mode_),
          ratio_mode_(other.ratio_mode_),
          old_(other.old_),
          pDistr_(other.pDistr_ ? other.pDistr_->clone_distrib() : nullptr),
          pProp_(other.pProp_ ? other.pProp_->clone() : nullptr)
    {
        if (pProp_)
            pProp_->rebind_distrib(pDistr_.get());
    }
    Variable::Variable(Variable &&) = default;

    Variable &Variable::operator=(Variable other)
    {
        swap(other);
        return *this;
    }

    void Variable::swap(Variable &other) noexcept
    {
        Attribute::swap(other);
        std::swap(pDistr_, other.pDistr_);
        std::swap(pProp_, other.pProp_);
    }
    void Variable::init_func_attr(std::string_view key, const Attribute &attr) { pFunc_->init_attr(key, attr); }

    void Variable::init_func_reshaping(bool discr)
    {
        if (has_model())
            discr = pDistr_->is_discrete();
        if (discr)
            pFunc_->init_discretization();
        pFunc_->init_truncation();

        // assert no missing arguments
        pFunc_->assert_is_callable_explicitly();
    }

    void Variable::init_distr(std::string_view opt)
    {
        if (na::is(opt))
            return;
        else
            pDistr_ = compose_distr_registries(opt, "distribution", {"continuos", "discrete"});
    }

    void Variable::init_distr_attr(std::string_view key, const Attribute &attr)
    {
        pDistr_->init_attr(key, attr);
        if (key == "loc" && pProp_)
            pProp_->associate_location(pDistr_->loc());
    }

    void Variable::init_distr_reshaping(bool is_discr)
    {
        if (!has_model())
            return;

        if (is_discr)
            pDistr_->init_discretization();
        pDistr_->init_truncation();

        pDistr_->assert_is_callable_explicitly();

        // check that non-constant variable has non-zero dimensional range
        helper_assert_non_degenerate_range_();
    }

    void Variable::helper_assert_non_degenerate_range_()
    {
        const Attribute &lower = pDistr_->lower();
        const Attribute &upper = pDistr_->upper();

        // skip if one of the range extremes is not provided
        if (!lower.is_callable() || !upper.is_callable())
            return;

        // skip if range is not degenerate
        const double sx = lower.value();
        const double dx = upper.value();
        if (sx < dx)
            return;

        // case empty interval, lower and upper are "open"
        if (sx > dx || !pDistr_->is_valid(sx) || !pDistr_->is_valid(dx))
            throw std::invalid_argument(
                "Definition produced an inconsistent range " + pDistr_->str_range() +
                " for " + err_name() + ".");

        if (!lower.is_constant() || !upper.is_constant())
            return;

        // warn user if it was not chosen as constant
        if (!is_user_const_)
        {
            msg::warning(
                "Collapsed (at least) one " + err_name() + " to constant.\n" +
                "Definitions of 'min' as " + lower.err_name() + " and 'max' as " + upper.err_name() +
                " produced a degenerate range.");
            init_numeric(sx);
            helper_init_null_proposal_();
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Variable : Setter (Proposal)  //////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void Variable::init_proposal(const std::string &proposal, double proposal_scale)
    {
        if (mode_ != Mode::Empty && mode_ != Mode::Numeric)
            return;

        if (na::is(proposal))
            return;

        if (mode_ == Mode::Empty)
            init_numeric(na::DOUBLE);
        is_user_const_ = false;

        if (proposal == "sample")
            helper_init_sample_proposal_();
        else
            helper_init_monte_carlo_proposal_(proposal, proposal_scale);
    }

    void Variable::helper_init_null_proposal_()
    {
        // assign the pointer to the proposer
        pProp_ = std::make_unique<NullProposal_>();

        // setup sample proposal functional methods
        ratio_mode_ = Ratio::None;
        valid_mode_ = Validity::Always;
    }

    void Variable::helper_init_sample_proposal_()
    {
        // assign the pointer to the proposer
        pProp_ = std::make_unique<SampleProposal_>(pDistr_.get());

        // setup sample proposal functional methods
        ratio_mode_ = Ratio::None;
        valid_mode_ = Validity::Always;
    }

    void Variable::helper_init_monte_carlo_proposal_(const std::string &proposal, double proposal_scale)
    {
        // raise a warning if "lnorm" proposal when variable can be negative
        if (pDistr_ && proposal == "lnorm" &&
            !pDistr_->lower().is_empty() && pDistr_->lower().value() < 0)
            msg::warning(
                "Provided a lognormal proposal for non-positive " + err_name() + ".");

        // select proposer
        pProp_ = proposal_switch(proposal, proposal_scale, is_discrete());

        // set modes
        ratio_mode_ = Ratio::Compute;
        valid_mode_ = pDistr_ ? Validity::Check : Validity::Always;
    }

    void Variable::user_reset(rnd_gen &prng)
    {
        if (is_constant())
            return;
        if (has_model())
            set(pDistr_->rand(prng));
        if (pProp_)
            pProp_->user_reset();
    }

    void Variable::optimize()
    {
        Attribute::optimize();
        if (pDistr_)
            pDistr_->optimize();
    }

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
        default:
            __builtin_unreachable();
        }
    }

    void Variable::clear() noexcept
    {
        mode_ = Mode::Empty;
        value_ = na::DOUBLE;
        pAttr_ = nullptr;
        pFunc_ = nullptr;
        pDistr_ = nullptr;
        pProp_ = nullptr;
    }

    bool Variable::is_discrete() const noexcept { return has_model() ? pDistr_->is_discrete() : Attribute::is_discrete(); }
    bool Variable::has_model() const noexcept { return pDistr_ != nullptr; }
    bool Variable::is_adaptable() const noexcept { return pProp_ && pProp_->is_adaptable(); }

    void Variable::propose(rnd_gen &prng) { set(pProp_->propose(prng, value_)); }
    bool Variable::is_valid_move() const
    {
        switch (valid_mode_)
        {
        case Validity::Always:
            return true;
        case Validity::Check:
            return pDistr_->is_valid(value_);
        default:
            __builtin_unreachable();
        }
    }

    double Variable::get_distr_log_pdf() const { return pDistr_->log_pdf(value_); }
    double Variable::get_distr_log_ratio() const { return pDistr_->log_ratio(old_, value_); }
    double Variable::get_propose_log_ratio() const { return pProp_->log_ratio(old_, value_); }

    void Variable::update_accept_rate(bool is_accepted) { pProp_->update_accept_rate(is_accepted); }
    void Variable::reset_accept_rate() { pProp_->reset_accept_rate(); }
    void Variable::adapt_proposal(double ratio) { pProp_->adapt(ratio); }

#ifndef NDEBUG
    std::string Variable::str_summary() const
    {
        return Attribute::str_summary() + (has_model() ? " | ~" + pDistr_->str_summary(false) : "");
    }
#endif
}
