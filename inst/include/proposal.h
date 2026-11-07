/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "tools/dictionary.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parser.h"
#include "tools/random.h"
#include "tools/string.h"
#include "distribution.h"
#include "function.h"
#include "registry.h"

namespace housepi
{
    //////////////////////////////  Proposal  //////////////////////////////

    /// Base of every Metropolis proposal.
    class Proposal
    {
    public:
        /// Default constructor.
        Proposal() = default;
        /// Virtual destructor.
        virtual ~Proposal() = default;

        /// Reset an adaptive proposal before moving, from `prior`.
        virtual void reset(const Distribution * /*prior*/) noexcept {}
        /// Adapt the proposal to acceptance probability `alpha` of the last move, with gain `gain`.
        virtual void adapt(double /*alpha*/, double /*gain*/) {}

        /// Propose a value from current value `xo`.
        virtual double propose(RndGen & /*prng*/, double /*xo*/) = 0;
        /// Compute the Metropolis-Hastings log-ratio of proposal densities, from `xo` to `xn`.
        virtual double log_ratio(double /*xo*/, double /*xn*/) const = 0;

        /// Name of the proposal as declared: `rw`, `rw_log`, `rw_unif`, `forward` or a distribution name.
        /// * the default and null proposals -> empty.
        virtual std::string name() const = 0;
        /// Coordinate the jumps are drawn in, empty when the kernel draws none.
        virtual std::string coord() const { return ""; }
        /// Current scale of the jumps, NA when the kernel has none.
        virtual double scale() const noexcept { return na::DOUBLE; }

        /// Whether propose() draws independently of the current value.
        /// Such a proposal can seed an observation-feasible initial value.
        virtual bool is_independent() const noexcept { return false; }

        /// Resolve the read-surface of an owned proposal distribution, evaluation only.
        /// @param owner Name of the variable proposed, for the error it may raise.
        virtual void finalize(const std::string & /*owner*/) {}
    };

    //////////////////////////////  NullProposal  //////////////////////////////

    /// Keep the current value: the variable never moves.
    class NullProposal final : public Proposal
    {
    public:
        NullProposal() {}

        std::string name() const override final { return ""; }

        double propose(RndGen & /*prng*/, double xo) override final { return xo; }
        double log_ratio(double /*xo*/, double /*xn*/) const override final { return 0.; }
    };

    //////////////////////////////  ForwardProposal  //////////////////////////////

    /// Draw from the prior itself.
    class ForwardProposal final : public Proposal
    {
    public:
        static constexpr const char *NAME = "forward";

        ForwardProposal(const Distribution *distr) : distr_(distr) {}

        std::string name() const override final { return NAME; }

        double propose(RndGen &prng, double /*xo*/) override final { return distr_->rand(prng); }
        double log_ratio(double xo, double xn) const override final { return distr_->log_ratio(xo, xn); }

    private:
        /// Prior drawn from.
        const Distribution *distr_;
    };

    //////////////////////////////  IndependenceProposal  //////////////////////////////

    /// Independence sampler drawing from an owned distribution, not the prior.
    class IndependenceProposal final : public Proposal
    {
    public:
        /// Build distribution `distr_name` with arguments `args`, truncated at the domain of `prior`.
        /// The domain is wired before reshaping.
        IndependenceProposal(
            const std::string &distr_name,
            const dict<Attribute> &args,
            const Distribution *prior)
            : distr_(registry::create_distribution(
                  distr_name,
                  ErrType::Proposal,
                  registry::distr_search_order(prior->is_discretized()))),
              prior_(prior)
        {
            for (const auto &[key, attr] : args)
                distr_->init_attr(key, attr);

            wire_to_prior_domain_();
            distr_->reshape(prior->is_discretized());
        }

        std::string name() const override final { return distr_->name(); }

        double propose(RndGen &prng, double /*xo*/) override final { return distr_->rand(prng); }
        double log_ratio(double xo, double xn) const override final { return distr_->log_ratio(xo, xn); }
        bool is_independent() const noexcept override final { return true; }

        /// Finalize the distribution, then assert it covers the domain of the prior on both sides.
        void finalize(const std::string &owner) override final
        {
            distr_->finalize();

            assert_covers_side_(intrinsic_min_, true, owner);
            assert_covers_side_(intrinsic_max_, false, owner);
        }

    private:
        /// Distribution drawn from.
        std::unique_ptr<Distribution> distr_;
        /// Prior, read for the domain to cover.
        const Distribution *prior_;

        /// Own bounds of the proposal, snapshot before it inherits the prior's.
        ///
        /// * the transfer runs at construction, ahead of reshape() settling the truncation mode;
        /// * it overwrites the bounds: finalize() checks coverage against this snapshot;
        /// * user arguments land first: a given `min` is captured here too.
        double intrinsic_min_ = na::DOUBLE;
        double intrinsic_max_ = na::DOUBLE;

        /// Snapshot the own support, then truncate the proposal at the domain of the prior.
        /// The domain may be live or data-bound: bounds stay live, so rand() re-truncates per move.
        void wire_to_prior_domain_()
        {
            intrinsic_min_ = distr_->min().is_callable() ? distr_->min().value() : na::DOUBLE;
            intrinsic_max_ = distr_->max().is_callable() ? distr_->max().value() : na::DOUBLE;

            // location shift of the prior: the proposal samples in the same absolute frame
            // bounds below stay in the shared delocated frame
            if (!prior_->loc().is_empty())
                distr_->init_attr("loc", prior_->loc());

            const Attribute &min_val = prior_->min();
            const Attribute &max_val = prior_->max();
            if (!min_val.is_empty())
                distr_->transfer_min(min_val);
            if (!max_val.is_empty())
                distr_->transfer_max(max_val);
        }

        /// Assert the own support of the proposal covers the domain of the prior on one side.
        /// The prior bound is read through its codomain: a live bound places every value it reaches.
        void assert_covers_side_(double intrinsic, bool below, const std::string &owner) const
        {
            // unbounded family: covers any domain on this side
            if (na::is(intrinsic))
                return;

            // an intrinsically bounded family covers only a prior bound it reaches in full
            // an unbounded prior, or a bound reaching past the family, leaves a gap
            const Codomain reach = (below ? prior_->min() : prior_->max()).codomain();
            if ((below ? reach.below : reach.above) == Interval::Unbounded ||
                (below ? intrinsic > reach.inf : intrinsic < reach.sup))
                raise_uncovered_(owner);
        }

        /// Throw: the support of the proposal cannot cover the domain of the prior.
        [[noreturn]] void raise_uncovered_(const std::string &owner) const
        {
            Codomain support;
            if (!na::is(intrinsic_min_))
                support.inf = intrinsic_min_, support.below = Interval::Closed;
            if (!na::is(intrinsic_max_))
                support.sup = intrinsic_max_, support.above = Interval::Closed;
            msg::stop("Proposal " + msg::val(distr_->name()) + " has support " + support.str() +
                      " that cannot cover the domain " + prior_->domain().str() +
                      " of " + owner + ".\nCorrect the proposal rule.");
        }
    };

    //////////////////////////////  Acceptance  //////////////////////////////

    /// Target acceptance probability of a walk of dimension `d`, its optimum on Gaussian targets.
    /// * `d` <= 10 -> acceptance maximizing the expected squared jump distance, from Theorem 1;
    /// * `d` > 10 -> 0.234, its limit as `d` grows.
    /// Based on @cite sherlock_optimal_2009
    constexpr double target_accept_rate(size_t d) noexcept
    {
        constexpr double RATES[] = {0.439, 0.351, 0.315, 0.296, 0.284, 0.276, 0.270, 0.266, 0.262, 0.259};
        return d > 10 ? 0.234 : RATES[d - 1];
    }

    /// Return acceptance probability min(1, e^`r`) of log-ratio `r`.
    /// * NaN -> 0.
    inline double accept_probability(double r)
    {
        return r >= 0. ? 1. : (r > -inf::DOUBLE ? std::exp(r) : 0.);
    }

    //////////////////////////////  AdaptableProposal  //////////////////////////////

    /// Base of random-walk proposals adapting their scale.
    class AdaptableProposal : public Proposal
    {
    public:
        //////  Constructor  //////

        /// Build continuous distribution `distr_name`, the law of the jumps.
        /// * discrete `prior` -> jumps rounded, see is_lattice();
        /// * NA `scale` -> adaptive;
        /// * otherwise -> fixed, an error unless positive and finite.
        AdaptableProposal(const std::string &distr_name, double scale, const Distribution *prior)
            : distr_(registry::create_distribution(distr_name, ErrType::Proposal, registry::distr_search_order(false))),
              is_lattice_(prior->is_discretized())
        {
            is_adaptive_ = na::is(scale);
            if (!is_adaptive_ && !(std::isfinite(scale) && scale > 0.))
                msg::stop("Proposal scale " + msg::val(str::format(scale)) + " is not a positive finite number.\nSet " +
                          msg::arg("proposal_scale") + " above 0, or leave it out to adapt the scale.");
            distr_->reshape(false);
        }

        /// Destructor.
        ~AdaptableProposal() override = default;

        //////  Adaptation  //////

        /// Reset before moving: default scale from `prior` if adaptive.
        void reset(const Distribution *prior) noexcept override final
        {
            if (is_adaptive())
                set_scale_(default_scale(prior));
        }

        /// Scale by Robbins-Monro multiplier e^(`gain` (`alpha` - target)), if adaptive.
        /// Based on @cite andrieu_tutorial_2008
        void adapt(double alpha, double gain) override final
        {
            if (is_adaptive_)
                set_scale_(get_scale_() * std::exp(gain * (alpha - TARGET_ACCEPT_RATE)));
        }

        /// Whether the scale is adaptive.
        bool is_adaptive() const noexcept { return is_adaptive_; }
        /// Whether the jumps are rounded half away from zero, onto the unit lattice of a discrete prior.
        /// Rounding is odd in the jump, so a symmetric jump law keeps a log-ratio of 0.
        bool is_lattice() const noexcept { return is_lattice_; }
        /// Current scale of the proposal distribution.
        double scale() const noexcept override final { return get_scale_(); }
        /// Return the default scale, from `prior`.
        virtual double default_scale(const Distribution * /*prior*/) const = 0;

        /// Resolve the read-surface of the owned proposal distribution.
        void finalize(const std::string & /*owner*/) override final { distr_->finalize(); }

    protected:
        static constexpr double TARGET_ACCEPT_RATE = target_accept_rate(1);

        std::unique_ptr<Distribution> distr_;
        /// Whether the scale is adaptive.
        bool is_adaptive_ = true;
        /// Whether the jumps are rounded, see is_lattice().
        const bool is_lattice_;

        /// Initialize scale at construction.
        void init_scale_(double scale) { set_scale_(is_adaptive() ? 1. : scale); }

        /// Read current scale.
        virtual double get_scale_() const = 0;
        /// Set scale of the proposal distribution.
        virtual void set_scale_(double /*scale*/) = 0;

        /// Return quartiles Q25 and Q75 of `prior`.
        static std::pair<double, double> get_iqr_from_prior_(const Distribution *prior)
        {
            return {prior->quantile(0.25), prior->quantile(0.75)};
        }
        /// Whether the interquartile range is finite and positive, fit for an initial scale.
        static bool is_iqr_finite_(double q25, double q75) noexcept
        {
            return std::isfinite(q25) && std::isfinite(q75) && q75 > q25;
        }
    };

    //////////////////////////////  RndWalkUnifProposal  //////////////////////////////

    /// Random walk with uniform jumps in [-scale, scale].
    /// * discrete prior -> jumps in [-scale - 1/2, scale + 1/2], rounded half away from zero.
    /// A whole `scale` jumps uniformly over the integers within it.
    class RndWalkUnifProposal final : public AdaptableProposal
    {
    public:
        static constexpr const char *NAME = "rw_unif";

        RndWalkUnifProposal(double scale, const Distribution *prior)
            : AdaptableProposal("uniform", scale, prior) { init_scale_(scale); }

        std::string name() const override final { return NAME; }
        std::string coord() const override final { return "identity"; }

        /// Jump `xo` by a uniform draw, rounded on a lattice.
        double propose(RndGen &prng, double xo) override final
        {
            const double dz = distr_->rand(prng);
            return xo + (is_lattice_ ? std::round(dz) : dz);
        }
        double log_ratio(double, double) const override final { return 0.; }

        /// Default jump half-width from `prior`: its interquartile range.
        /// * infinite or empty range -> 1.
        double default_scale(const Distribution *prior) const override final
        {
            const auto [q25, q75] = get_iqr_from_prior_(prior);
            if (is_iqr_finite_(q25, q75))
                return q75 - q25;
            return 1.;
        }

    protected:
        double get_scale_() const override final { return distr_->max().value() - (is_lattice_ ? 0.5 : 0.); }
        /// Set jumps to [-scale, scale], widened by 1/2 on a lattice.
        void set_scale_(double scale) override final
        {
            const double s = scale + (is_lattice_ ? 0.5 : 0.);
            distr_->set_min(-s);
            distr_->set_max(s);
        }
    };

    //////////////////////////////  RndWalkProposal  //////////////////////////////

    /// Random walk with normal jumps in a coordinate z of the value, scale as standard deviation.
    /// * `Identity` -> z = x;
    /// * `Log` -> z = log(x - loc), `loc` the location of the prior;
    /// * `Half` -> z = log|x - c|, `c` the bound of a half-line domain;
    /// * `Logit` -> z = logit((x - a) / (b - a)), on a domain bounded by `a` < `b`;
    /// * `Lattice` -> z = x on a discrete prior, the jump rounded half away from zero;
    /// * `CheckedLattice` -> as `Lattice`, for the default walk.
    /// The jumps are drawn continuous in z.
    /// Bounds are read live: `Half` and `Logit` stay in the domain, `CheckedLattice` refuses a value outside it.
    class RndWalkProposal final : public AdaptableProposal
    {
    public:
        static constexpr const char *NAME = "rw";
        static constexpr const char *LOG_NAME = "rw_log";

        /// Coordinate the jumps are drawn in.
        enum class Coord : std::uint8_t
        {
            Identity,
            Log,
            Half,
            Logit,
            CheckedLattice,
            Lattice
        };

        /// Build declared walk `rw`, in `Identity`, or `rw_log`, in `Log`.
        /// * `rw` on a discrete prior -> `Lattice`;
        /// * NA `scale` -> adaptive.
        RndWalkProposal(double scale, const Distribution *prior, Coord coord)
            : AdaptableProposal("normal", scale, prior),
              coord_(coord == Coord::Identity && is_lattice_ ? Coord::Lattice : coord),
              is_declared_(true),
              domain_(&prior->domain())
        {
            if (coord == Coord::Log)
            {
                if (is_lattice_)
                    msg::stop("Proposal " + msg::val(LOG_NAME) + " cannot move on a discrete prior.\nUse " +
                              msg::val(NAME) + " or " + msg::val(RndWalkUnifProposal::NAME) + " instead.");

                // jumps stay above the location: the domain must too
                if (prior->min().is_empty() || prior->min().value() < 0)
                    msg::stop("Proposal " + msg::val(LOG_NAME) + " needs a prior with a non-negative minimum.\nSet " +
                              msg::arg("min") + " to 0 or above, or use " + msg::val(NAME) + " instead.");
            }
            init_scale_(scale);
        }

        /// Build the default walk, adaptive, in the coordinate mapping the domain of `prior` onto the real line.
        /// * discrete prior -> `CheckedLattice`.
        explicit RndWalkProposal(const Distribution *prior)
            : AdaptableProposal("normal", na::DOUBLE, prior),
              domain_(&prior->domain())
        {
            const bool lower = domain_->below != Interval::Unbounded;
            const bool upper = domain_->above != Interval::Unbounded;
            coord_ = is_lattice_ ? Coord::CheckedLattice
                                 : (lower && upper ? Coord::Logit : (lower || upper ? Coord::Half : Coord::Identity));
            side_ = upper && !lower ? -1. : 1.;
            init_scale_(na::DOUBLE);
        }

        std::string name() const override final
        {
            if (!is_declared_)
                return "";
            if (coord_ == Coord::Log)
                return LOG_NAME;
            return NAME;
        }
        std::string coord() const override final
        {
            switch (coord_)
            {
            case Coord::Log:
            case Coord::Half:
                return "log";
            case Coord::Logit:
                return "logit";
            // cases Coord::Identity, Coord::Lattice, Coord::CheckedLattice
            default:
                return "identity";
            }
        }

        /// Jump `xo` by a normal draw in its coordinate (see jump()).
        double propose(RndGen &prng, double xo) override final { return jump(xo, distr_->rand(prng)); }

        /// Jump `xo` by `dz` in its coordinate, and cache the log-ratio of the move.
        /// * log coordinate -> J(xn) - J(xo) = -dz, with J = log|dz/dx|;
        /// * lattice -> `dz` rounded, the log-ratio 0;
        /// * value rounded onto a bound, or `CheckedLattice` value outside the domain -> infinity, the move refused.
        double jump(double xo, double dz)
        {
            switch (coord_)
            {
            case Coord::Log:
                return log_step_(xo, dz, domain_->relocate(0.));
            case Coord::Half:
                return log_step_(xo, dz, bound_());
            case Coord::Logit:
            {
                // e^-z read from x, no logarithm of the current value
                const double a = lower_(), b = upper_();
                const double eo = (b - xo) / (xo - a), en = eo * std::exp(-dz);
                const double xn = a + (b - a) / (1. + en);
                log_ratio_ = xn > a && xn < b ? dz + 2. * std::log((1. + en) / (1. + eo)) : inf::DOUBLE;
                return xn;
            }
            case Coord::Lattice:
                log_ratio_ = 0.;
                return xo + std::round(dz);
            case Coord::CheckedLattice:
            {
                const double xn = xo + std::round(dz);
                log_ratio_ = domain_->is_valid(xn) ? 0. : inf::DOUBLE;
                return xn;
            }
            // case Coord::Identity:
            default:
                log_ratio_ = 0.;
                return xo + dz;
            }
        }
        /// Return the log-ratio cached by jump().
        /// @note Read after jump(), on the move it proposed.
        double log_ratio(double /*xo*/, double /*xn*/) const override final { return log_ratio_; }

        /// Map `x` to its coordinate.
        double to_coord(double x) const
        {
            switch (coord_)
            {
            case Coord::Log:
                return std::log(x - domain_->relocate(0.));
            case Coord::Half:
                return std::log(side_ * (x - bound_()));
            case Coord::Logit:
                return std::log((x - lower_()) / (upper_() - x));
            // cases Coord::Identity, Coord::Lattice, Coord::CheckedLattice
            default:
                return x;
            }
        }

        /// Default jump scale from `prior`: its interquartile range in the coordinate, over 1.349.
        /// * infinite or empty range -> 1.
        double default_scale(const Distribution *prior) const override final
        {
            const auto [q25, q75] = get_iqr_from_prior_(prior);
            // a coordinate decreasing in x swaps the quartiles
            double z25 = to_coord(q25), z75 = to_coord(q75);
            if (z25 > z75)
                std::swap(z25, z75);
            if (is_iqr_finite_(z25, z75))
                return (z75 - z25) / 1.349;
            return 1.;
        }

    protected:
        /// Coordinate the jumps are drawn in.
        Coord coord_ = Coord::Identity;
        /// Whether declared as `rw` or `rw_log`, false for the default walk.
        bool is_declared_ = false;
        /// Domain of the prior.
        const Domain *domain_;
        /// Side of the domain from its anchor: 1 above, -1 below.
        double side_ = 1.;
        /// Log-ratio of the last move proposed.
        double log_ratio_ = 0.;

        /// Lower bound, in the absolute frame.
        double lower_() const { return domain_->relocate(domain_->min.value()); }
        /// Upper bound, in the absolute frame.
        double upper_() const { return domain_->relocate(domain_->max.value()); }
        /// Bound of a half-line domain, in the absolute frame.
        double bound_() const { return side_ > 0. ? lower_() : upper_(); }

        /// Jump `xo` by `dz` in coordinate log|x - `c`|, and cache the log-ratio of the move.
        double log_step_(double xo, double dz, double c)
        {
            const double xn = c + (xo - c) * std::exp(dz);
            log_ratio_ = side_ * (xn - c) > 0. ? -dz : inf::DOUBLE;
            return xn;
        }

        double get_scale_() const override final { return distr_->arg(1); }
        /// Set standard deviation of the jumps.
        void set_scale_(double scale) override final { distr_->set_arg(1, scale); }
    };

    //////////////////////////////  Declared proposal  //////////////////////////////

    /// Whether proposal `name` is a random walk, tuned by a `scale` argument.
    inline bool is_random_walk(std::string_view name) noexcept
    {
        return name == RndWalkProposal::NAME ||
               name == RndWalkProposal::LOG_NAME ||
               name == RndWalkUnifProposal::NAME;
    }

    /// Names a declared proposal takes besides a distribution name.
    inline std::vector<std::string> kernel_names()
    {
        return {ForwardProposal::NAME, RndWalkProposal::NAME, RndWalkProposal::LOG_NAME, RndWalkUnifProposal::NAME};
    }

    /// Whether declared proposal `prop` may adapt its scale on a variable it builds.
    /// * none declared -> adaptive;
    /// * random walk, no `scale` -> adaptive;
    /// * random walk, numeric `scale` -> frozen;
    /// * random walk, `scale` resolved per variable -> adaptive where it resolves to NA;
    /// * any other proposal -> frozen.
    inline bool is_declared_adaptive(const parse::Definition &prop)
    {
        if (prop.is_na())
            return true;
        if (!is_random_walk(prop.root.func))
            return false;
        for (const parse::Entry &child : prop.root.children)
            if (child.key == "scale")
                return !(child.is_terminal() && parse::is_numeric(child.type));
        return true;
    }

    //////////////////////////////  block_rw  //////////////////////////////

    /// Adaptive Metropolis: block random walk with jumps from the Cholesky factor of the running covariance.
    /// Robbins-Monro tunes its global scale.
    /// Based on @cite haario_adaptive_2001 and @cite andrieu_tutorial_2008
    namespace block_rw
    {
        /// Fraction of the adaptation moving members alone.
        constexpr double WARMUP = 0.5;
        /// Steps of the first moment window after the warm-up, each next doubling.
        constexpr std::uint64_t WINDOW = 50;

        /// log(2.38), the log of the optimal random-walk scale 2.38 / sqrt(d) at d = 1.
        /// Based on @cite roberts_weak_1997
        constexpr double LOG_2_38 = 0.8671004876833833;
        /// Variance of rounding to a unit lattice.
        constexpr double LATTICE_VAR = 1. / 12.;

        /// Index of entry (`i`, `j`), `j` <= `i`, in a lower triangle packed by rows.
        constexpr size_t packed(size_t i, size_t j) noexcept { return i * (i + 1) / 2 + j; }

        /// Replace packed symmetric `a` with its lower Cholesky factor.
        /// @return Whether `a` is positive definite, `a` left partly overwritten otherwise.
        inline bool cholesky(std::vector<double> &a, size_t d)
        {
            for (size_t j = 0; j < d; ++j)
            {
                double s = a[packed(j, j)];
                for (size_t k = 0; k < j; ++k)
                    s -= a[packed(j, k)] * a[packed(j, k)];
                if (!(s > 0.))
                    return false;
                const double ljj = std::sqrt(s);
                a[packed(j, j)] = ljj;
                for (size_t i = j + 1; i < d; ++i)
                {
                    double t = a[packed(i, j)];
                    for (size_t k = 0; k < j; ++k)
                        t -= a[packed(i, k)] * a[packed(j, k)];
                    a[packed(i, j)] = t / ljj;
                }
            }
            return true;
        }

        /// Factor of one block, its running moments and the last proposal it drew.
        struct Factor
        {
            /// Lower factor L, packed by rows.
            std::vector<double> lower = {};
            /// Jump e^log_scale L u of the last proposal.
            std::vector<double> jump = {};
            /// Coordinates of the current state, filled before observe().
            std::vector<double> z = {};
            /// Running mean of the coordinates.
            std::vector<double> mean = {};
            /// Running sums of centred products of the coordinates, packed.
            std::vector<double> m2 = {};
            /// Scratch of the draw and of the moments.
            std::vector<double> work = {};
            /// Scratch of the factorization.
            std::vector<double> chol = {};
            /// Whether each member moves on a unit lattice.
            std::vector<bool> lattice = {};
            /// States observed.
            std::uint64_t n_obs = 0;
            /// Log of the multiplier of L.
            double log_scale = 0.;
            /// Whether every member moves on a unit lattice.
            bool all_lattice = false;

            /// Dimension of the block.
            size_t dim() const noexcept { return jump.size(); }
            /// Whether the last jump rounds to 0 on every member, all on a unit lattice.
            /// * some member off a lattice -> false, the jump never read.
            bool is_null_jump() const
            {
                if (!all_lattice)
                    return false;
                for (double dz : jump)
                    if (!(std::fabs(dz) < 0.5))
                        return false;
                return true;
            }
            /// Covariance e^(2 log_scale) L L^T of the joint jump, entry (`i`, `j`).
            double jump_cov(size_t i, size_t j) const
            {
                double acc = 0.;
                for (size_t k = 0; k <= std::min(i, j); ++k)
                    acc += lower[packed(i, k)] * lower[packed(j, k)];
                return std::exp(2. * log_scale) * acc;
            }
            /// Standard deviation of the jump of member `i`.
            double jump_sd(size_t i) const { return std::sqrt(jump_cov(i, i)); }

            /// Set L diagonal from `scales`, the multiplier to 2.38 / sqrt(d), and restart the adaptation.
            /// @param on_lattice Whether each member moves on a unit lattice.
            void init(const std::vector<double> &scales, std::vector<bool> on_lattice)
            {
                const size_t d = scales.size(), n_packed = d * (d + 1) / 2;
                lattice = std::move(on_lattice);
                all_lattice = std::find(lattice.begin(), lattice.end(), false) == lattice.end();
                lower.assign(n_packed, 0.);
                for (size_t i = 0; i < d; ++i)
                    lower[packed(i, i)] = scales[i];
                jump.assign(d, 0.);
                z.assign(d, 0.);
                mean.assign(d, 0.);
                m2.assign(n_packed, 0.);
                work.assign(d, 0.);
                chol.assign(n_packed, 0.);
                n_obs = 0;
                log_scale = LOG_2_38 - 0.5 * std::log(d);
            }

            /// Draw u ~ N(0, I) and set the jump e^log_scale L u.
            void propose(RndGen &prng)
            {
                const size_t d = dim();
                const double mult = std::exp(log_scale);
                for (size_t i = 0; i < d; ++i)
                    work[i] = standard_normal(prng);
                for (size_t i = 0; i < d; ++i)
                {
                    double acc = 0.;
                    for (size_t j = 0; j <= i; ++j)
                        acc += lower[packed(i, j)] * work[j];
                    jump[i] = mult * acc;
                }
            }

            /// Add coordinates `z` to the running moments.
            void observe()
            {
                const size_t d = dim();
                ++n_obs;
                for (size_t i = 0; i < d; ++i)
                {
                    work[i] = z[i] - mean[i];
                    mean[i] += work[i] / n_obs;
                }
                for (size_t i = 0; i < d; ++i)
                    for (size_t j = 0; j <= i; ++j)
                        m2[packed(i, j)] += work[i] * (z[j] - mean[j]);
            }

            /// Drop the running moments.
            void forget()
            {
                mean.assign(mean.size(), 0.);
                m2.assign(m2.size(), 0.);
                n_obs = 0;
            }

            /// Set L to the Cholesky factor of the running covariance, regularized.
            /// * lattice member -> covariance of its continuous embedding, LATTICE_VAR added;
            /// * fewer than 5 d states, or not positive definite -> L kept.
            void reshape()
            {
                const size_t d = dim();
                if (n_obs < 5 * d)
                    return;
                const double n1 = n_obs - 1;
                double tr = 0.;
                for (size_t i = 0; i < d; ++i)
                    tr += m2[packed(i, i)];
                const double eps = 1e-6 * tr / (d * n1) + 1e-12;
                for (size_t i = 0; i < d; ++i)
                    for (size_t j = 0; j <= i; ++j)
                        chol[packed(i, j)] = m2[packed(i, j)] / n1 + (i == j ? eps + lattice[i] * LATTICE_VAR : 0.);
                if (cholesky(chol, d))
                    lower.swap(chol);
            }

            /// Move log_scale by `gain` (`alpha` - target_accept_rate(d)).
            /// @param alpha Acceptance probability of the last proposal.
            void adapt(double alpha, double gain) { log_scale += gain * (alpha - target_accept_rate(dim())); }
        };
    }
}
