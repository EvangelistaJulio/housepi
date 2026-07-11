/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>

#include "tools/messages.h"
#include "tools/random.h"

#include "distribution.h"
#include "function.h"

namespace housepi
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

        /// If adaptive, reset to default proposal scale.
        virtual void reset_scale(const Distribution * /*prior*/) {}

        /// Check if proposal is adaptive.
        virtual bool is_adaptive() const noexcept { return false; }

        /// Reset the acceptance rate estimate.
        virtual void reset_accept_rate() noexcept {}
        /// Adapt parametric proposal based on current steps ratio.
        virtual void adapt(double /*ratio*/) {}
        /// Adapt parametric proposal based on current steps ratio.
        virtual void recover() noexcept {}

        /// Propose a value using random number generator.
        virtual double propose(rnd_gen & /*prng*/, double /*xo*/) = 0;
        /// Compute log-ratio of proposal densities for Metropolis-Hastings.
        virtual double log_ratio(double /*xo*/, double /*xn*/) const = 0;

        /// Rebind proposal to a new prior distribution.
        virtual void rebind_prior(const Distribution * /*new_distr*/) {}

        /// Resolve the read-surface of an owned proposal distribution (evaluation-only).
        virtual void optimize() {}
    };

#define MAKE_PROPOSAL_CLONE_METHOD(Class_) \
    std::unique_ptr<Proposal> clone() const override { return std::make_unique<Class_>(*this); }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  *** NullProposal : class ***  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class NullProposal final : public Proposal
    {
    public:
        NullProposal() {}
        MAKE_PROPOSAL_CLONE_METHOD(NullProposal)

        double propose(rnd_gen & /*prng*/, double xo) override final { return xo; }
        double log_ratio(double /*xo*/, double /*xn*/) const override final { return 0.; }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** PriorProposal : class ***  /////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class PriorProposal final : public Proposal
    {
    public:
        PriorProposal(const Distribution *distr) : pDistr_(distr) {}
        MAKE_PROPOSAL_CLONE_METHOD(PriorProposal)

        void rebind_prior(const Distribution *new_distr) override final { pDistr_ = new_distr; }

        double propose(rnd_gen &prng, double /*xo*/) override final { return pDistr_->rand(prng); }
        double log_ratio(double xo, double xn) const override final { return pDistr_->log_ratio(xo, xn); }

    private:
        const Distribution *pDistr_;
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** IndependenceProposal : class ***  /////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// Independence sampler drawing from an owned distribution (not the prior).
    class IndependenceProposal final : public Proposal
    {
    public:
        IndependenceProposal(std::unique_ptr<Distribution> distr) : pDistr_(std::move(distr)) {}
        IndependenceProposal(const IndependenceProposal &other)
            : pDistr_(other.pDistr_ ? other.pDistr_->clone_distrib() : nullptr) {}
        MAKE_PROPOSAL_CLONE_METHOD(IndependenceProposal)

        double propose(rnd_gen &prng, double /*xo*/) override final { return pDistr_->rand(prng); }
        double log_ratio(double xo, double xn) const override final { return pDistr_->log_ratio(xo, xn); }

        void optimize() override final
        {
            if (pDistr_)
                pDistr_->optimize();
        }

    private:
        std::unique_ptr<Distribution> pDistr_; /// distribution drawn from independently
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

        AdaptableProposal(
            const std::string &distr_name,
            double scale,
            bool is_discrete)
            : pDistr_(compose_distr_registries(distr_name, "proposal", {"continuos"})),
              is_adaptive_(na::is(scale))
        {
            if (is_discrete)
                pDistr_->init_discretization();
            pDistr_->init_truncation();
        }
        virtual ~AdaptableProposal() = default;

        AdaptableProposal(const AdaptableProposal &other)
            : pDistr_(other.pDistr_ ? other.pDistr_->clone_distrib() : nullptr),
              is_adaptive_(other.is_adaptive_) {}

        //////////////////////////
        //////  ADAPTATION  //////
        //////////////////////////

        /// NOTE: initialised at the target to avoid first steps adaptation.
        void reset_accept_rate() noexcept override final { accept_rate_ = TARGET_ACCEPT_RATE; }

        void reset_scale(const Distribution *prior) override final
        {
            if (is_adaptive())
                set_scale_(get_default_scale_(prior));
        }

        void adapt(double ratio) override final
        {
            set_scale_(get_scale_() * robbins_monro_(accept_rate_, ratio));
        }

        /// Wheter proposal is adaptive.
        bool is_adaptive() const noexcept override final { return is_adaptive_; }

        /// Resolve the read-surface of the owned proposal distribution.
        void optimize() override final
        {
            if (pDistr_)
                pDistr_->optimize();
        }

        /// Correct acceptance rate for refused move.
        /// NOTE: exponential moving average for accept = FALSE.
        void recover() noexcept override final { accept_rate_ -= EMA_SMOOTHING; }

        /// Set EMA smoothing factor (0, 1). NOTE: shared across variables.
        static void set_smoothing(double s) noexcept { EMA_SMOOTHING = s; }
        /// Get EMA smoothing factor.
        static double get_smoothing() noexcept { return EMA_SMOOTHING; }

    protected:
        /// Target Metropolis acceptance rate.
        /// Based on @cite roberts_weak_1997
        static constexpr double TARGET_ACCEPT_RATE = 0.44;
        /// EMA smoothing factor. NOTE: effective window is 1. / EMA_SMOOTHING.
        inline static double EMA_SMOOTHING = 1.;

        std::unique_ptr<Distribution> pDistr_;    /// distribution associated to proposal
        bool is_adaptive_ = true;                 /// whether scale is adaptive
        double accept_rate_ = TARGET_ACCEPT_RATE; /// acceptance rate

        /// Set default proposal scale.
        virtual double get_default_scale_(const Distribution * /*prior*/) = 0;
        /// Read the current scale from pDistr_.
        virtual double get_scale_() const = 0;
        /// Mechanism to change corresponding scale in pDistr_.
        virtual void set_scale_(double /*scale*/) = 0;

        /// Exponentially moving average (EMA) of the acceptance indicator.
        /// NOTE: exponential moving average for accept = TRUE.
        void forward_accept_rate_() { accept_rate_ += EMA_SMOOTHING * (1. - accept_rate_); }

        /// Vanishing Robbins-Monro multiplier from an observed accept rate.
        /// Based on @cite andrieu_tutorial_2008
        static double robbins_monro_(double accept_rate, double ratio) noexcept
        {
            return std::exp((1. - ratio) * (accept_rate - TARGET_ACCEPT_RATE));
        }
        /// Compute InterQuartile Range (IQR) from prior distribution (Q25, Q75).
        static std::pair<double, double> get_iqr_from_prior_(const Distribution *prior)
        {
            return {prior->quantile(0.25), prior->quantile(0.75)};
        }
        /// Wheter IQR is finite and can be used to compute initial proposal scale.
        static bool is_iqr_finite_(double q25, double q75) noexcept
        {
            return std::isfinite(q25) && std::isfinite(q75) && q75 > q25;
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** RndWalkUnifProposal : class ***  //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class RndWalkUnifProposal final : public AdaptableProposal
    {
    public:
        RndWalkUnifProposal(double scale, const Distribution *prior)
            : AdaptableProposal("uniform", scale, prior->is_discrete())
        {
            set_scale_(is_adaptive() ? 1. : std::fabs(scale));
        }
        MAKE_PROPOSAL_CLONE_METHOD(RndWalkUnifProposal)

        double propose(rnd_gen &prng, double xo) override final
        {
            forward_accept_rate_();
            return xo + pDistr_->rand(prng);
        }
        double log_ratio(double, double) const override final { return 0.; }

    protected:
        double get_default_scale_(const Distribution *prior) override final
        {
            const auto [q25, q75] = get_iqr_from_prior_(prior);
            if (is_iqr_finite_(q25, q75))
                return q75 - q25;
            return 1.;
        }
        /// NOTE: uniform stores its scale in the bounds, not in the args.
        double get_scale_() const override final { return pDistr_->upper().value(); }
        /// Adapt standard minimum and maximum.
        void set_scale_(double scale) override final
        {
            pDistr_->set_lower(-scale);
            pDistr_->set_upper(scale);
        }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  *** RndWalkProposal : class ***  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class RndWalkProposal final : public AdaptableProposal
    {
    public:
        RndWalkProposal(double scale, const Distribution *prior)
            : AdaptableProposal("normal", scale, prior->is_discrete())
        {
            set_scale_(is_adaptive() ? 1. : std::fabs(scale));
        }
        MAKE_PROPOSAL_CLONE_METHOD(RndWalkProposal)

        double propose(rnd_gen &prng, double xo) override final
        {
            forward_accept_rate_();
            return xo + pDistr_->rand(prng);
        }
        double log_ratio(double, double) const override final { return 0.; }

    protected:
        double get_default_scale_(const Distribution *prior) override final
        {
            const auto [q25, q75] = get_iqr_from_prior_(prior);
            if (is_iqr_finite_(q25, q75))
                return (q75 - q25) / 1.349;
            return 1.;
        }
        double get_scale_() const override final { return pDistr_->arg(1); }
        /// Adapt standard deviation.
        void set_scale_(double scale) override final { pDistr_->set_arg(1, scale); }
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  *** RndWalkLogProposal : class ***  //////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    class RndWalkLogProposal final : public AdaptableProposal
    {
    public:
        RndWalkLogProposal(double scale, const Distribution *prior)
            : AdaptableProposal("lognormal", scale, prior->is_discrete())
        {
            // assert that value range is positive
            if (prior->lower().is_empty() || prior->lower().value() < 0)
                msg::stop(
                    "Provided proposal " + msg::val("lognormal") + " for non-positive");

            // associate location parameter from prior
            if (prior->loc().is_callable())
                pAttr_loc_ = prior->loc().pointer();
            mode_ = pAttr_loc_ ? Location::Active : Location::None;

            set_scale_(is_adaptive() ? 1. : std::fabs(scale));
        }
        MAKE_PROPOSAL_CLONE_METHOD(RndWalkLogProposal)

        double propose(rnd_gen &prng, double xo) override final
        {
            forward_accept_rate_();
            switch (mode_)
            {
            case Location::None:
                return xo * pDistr_->rand(prng);
            default: // Location::Active:
                const double loc = pAttr_loc_->value();
                return loc + (xo - loc) * pDistr_->rand(prng);
            }
        }
        double log_ratio(double xo, double xn) const override final
        {
            switch (mode_)
            {
            case Location::None:
                return std::log(xo / xn);
            default: // Location::Active:
                const double loc = pAttr_loc_->value();
                return std::log((xo - loc) / (xn - loc));
            }
        }

    protected:
        enum class Location : std::uint8_t
        {
            None,
            Active
        } mode_ = Location::None;
        const Attribute *pAttr_loc_ = nullptr;

        double get_default_scale_(const Distribution *prior) override final
        {
            const double loc = pAttr_loc_ ? pAttr_loc_->value() : 0.;
            const auto [q25, q75] = get_iqr_from_prior_(prior);
            if (is_iqr_finite_(q25 - loc, q75 - loc))
                return std::log((q75 - loc) / (q25 - loc)) / 1.349;
            return 1.;
        }
        double get_scale_() const override final { return pDistr_->arg(1); }
        /// Adapt log of standard deviation.
        void set_scale_(double scale) override final { pDistr_->set_arg(1, scale); }
    };
}