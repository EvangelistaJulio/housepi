/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/gamma.hpp>

/// Block boost long double policy.
using no_promote_policy = boost::math::policies::policy<
    boost::math::policies::promote_float<false>,
    boost::math::policies::promote_double<false>>;

/// Discrete quantile policy of R: min{k : P(X <= k) >= p}.
///
/// * boost's default `integer_round_outwards` rounds toward the tails, unlike R for p < 0.5;
/// * `integer_round_up` matches R exactly.
using discrete_up_policy = boost::math::policies::policy<
    boost::math::policies::promote_float<false>,
    boost::math::policies::promote_double<false>,
    boost::math::policies::discrete_quantile<boost::math::policies::integer_round_up>>;

#include "tools/math.h"
#include "tools/messages.h"
#include "tools/string.h"

#include "distribution.h"
#include "registry.h"

namespace housepi
{
    //////////////////////////////  registry::RegDistributions  //////////////////////////////

    namespace registry
    {
        RegDistributions::func_uptr RegDistributions::create(std::string_view name,
                                                             bool is_cumul) const
        {
            auto it = factory_map_.find(name);
            if (it == factory_map_.end())
                return nullptr;

            auto output = it->second();
            if (is_cumul)
                output->set_eval_as_cdf();
            else
                output->set_eval_as_pdf();

            return output;
        }

        std::unique_ptr<Distribution> create_distribution(std::string_view name,
                                                          ErrType err_type,
                                                          const SearchOrder &order,
                                                          bool is_cumul)
        {
            if (na::is(name))
                return nullptr;
            return std::unique_ptr<Distribution>(
                static_cast<Distribution *>(create_function(name, err_type, order, is_cumul).release()));
        }
    }
}

//////////////////////////////  Distribution implementations  //////////////////////////////

namespace housepi
{
    namespace
    {
        using std::exp;
        using std::lgamma;
        using std::log;
        using std::log1p;
        using std::pow;
        using std::sqrt;

        /// Logarithm of binomial coefficient.
        double log_binom_coeff(unsigned N, unsigned k)
        {
            if (N == 0 || k == N)
                return 0.;
            return lgamma(N + 1) - lgamma(k + 1) - lgamma(N - k + 1);
        }

        /// Log-ratio of binomial coefficients.
        double logratio_binom_coeff(unsigned N, unsigned ko, unsigned kn)
        {
            if (N == 0 || ko == kn)
                return 0;
            return lgamma(ko + 1) + lgamma(N - ko + 1) - lgamma(kn + 1) - lgamma(N - kn + 1);
        }
        /// Logarithm of beta function.
        double lbeta(double a, double b)
        {
            return lgamma(a) + lgamma(b) - lgamma(a + b);
        }

        /// Value derived from arguments, stored by finalize() when the arguments it reads are constant.
        /// * stored -> read per evaluation;
        /// * NA -> computed per evaluation.
        struct ArgCache
        {
            double value = na::DOUBLE; /// Stored value, NA when not stored.

            /// Whether the value is stored.
            bool is_stored() const noexcept { return !na::is(value); }
            /// Store `compute()` when `is_constant`, NA otherwise.
            template <class F>
            void store(bool is_constant, F compute) { value = is_constant ? compute() : na::DOUBLE; }
            /// Stored value, `compute()` when not stored.
            template <class F>
            double read(F compute) const { return is_stored() ? value : compute(); }
        };

        /// Density at the origin of a law behaving as x^(shape - 1) near it.
        /// * shape below 1 -> infinity;
        /// * shape 1 -> 1 / scale;
        /// * shape above 1 -> 0.
        double origin_pdf(double shape, double scale)
        {
            if (shape == 1.)
                return 1. / scale;
            return shape < 1. ? inf::DOUBLE : 0.;
        }

        /// Define left discretization, with a minimum closed at 0.
#define INIT_DISCR_LEFT_MIN_ZERO_CLOSE()        \
    void discretize() override final            \
    {                                           \
        discr_mode_ = Discrete::Left;           \
    }                                           \
    void truncate() override final              \
    {                                           \
        if (helper_init_domain_min_(0., false)) \
            helper_init_distr_truncation_();    \
    }

        /// Define a discrete-defined distribution, with a minimum closed at 0.
#define INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()     \
    void discretize() override final {}         \
    void truncate() override final              \
    {                                           \
        if (helper_init_domain_min_(0., false)) \
            helper_init_distr_truncation_();    \
    }

        /// Refuse discretization, for a law with no discrete counterpart.
#define INIT_DISCR_REFUSED()                                             \
    void discretize() override final                                     \
    {                                                                    \
        msg::stop("Invoked invalid option discrete for " + str() + "."); \
    }

        //////////////////////////////  Uniform  //////////////////////////////

        class Uniform_ final : public Cloneable<Uniform_, Distribution>
        {
        public:
            static constexpr const char *NAME = "uniform";
            Uniform_() : Cloneable({}, 0., 1.) {}

        private:
            double norm_() const { return domain_.max.value() - domain_.min.value(); }

        protected:
            double raw_function_(double) const override final { return 1. / norm_(); }
            double raw_primitive_(double x) const override final { return (x - domain_.min.value()) / norm_(); }
            double raw_ccdf_(double x) const override final { return (domain_.max.value() - x) / norm_(); }
            double raw_log_func_(double) const override final { return -log(norm_()); }
            double raw_log_ratio_(double, double) const override final { return 0.; }
            double raw_rand_(RndGen &prng) const override final
            {
                const double r = standard_uniform(prng);
                return domain_.min.value() * (1. - r) + domain_.max.value() * r;
            }
            double raw_quantile_(double p) const override final
            {
                return domain_.min.value() * (1. - p) + domain_.max.value() * p;
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Uniform_, Continuous)

        //////////////////////////////  Loguniform  //////////////////////////////

        /// Reciprocal distribution: flat in log x over [min, max], 0 < min < max.
        /// * support -> the domain itself, no truncation;
        /// * `min` and `max` -> required, no default.
        class Loguniform_ final : public Cloneable<Loguniform_, Distribution>
        {
        public:
            static constexpr const char *NAME = "loguniform";
            Loguniform_() = default;

            INIT_DISCR_REFUSED()

        private:
            /// Log-width log(max / min).
            double norm_() const { return log1p((domain_.max.value() - domain_.min.value()) / domain_.min.value()); }

        protected:
            /// Both bounds positive, so both required.
            ArgRequirements arg_requirements_() const override final
            {
                return {{&domain_.min, Codomain::positive()}, {&domain_.max, Codomain::positive()}};
            }

            double raw_function_(double x) const override final { return 1. / (x * norm_()); }
            double raw_primitive_(double x) const override final
            {
                const double a = domain_.min.value();
                return log1p((x - a) / a) / norm_();
            }
            double raw_ccdf_(double x) const override final { return log1p((domain_.max.value() - x) / x) / norm_(); }
            double raw_log_func_(double x) const override final { return -log(x * norm_()); }
            double raw_log_ratio_(double xo, double xn) const override final { return log(xo / xn); }
            double raw_quantile_(double p) const override final
            {
                const double a = domain_.min.value();
                return a + a * math::expm1(p * norm_());
            }
            double raw_rand_(RndGen &prng) const override final { return raw_quantile_(standard_uniform(prng)); }
        };
        HOUSEPI_REGISTER_FUNCTION(Loguniform_, Continuous)

        //////////////////////////////  Discrete uniform  //////////////////////////////

        /// Flat mass over the integers in [min, max].
        class UniformDiscrete_ final : public Cloneable<UniformDiscrete_, Distribution>
        {
        public:
            static constexpr const char *NAME = "uniform";
            UniformDiscrete_() : Cloneable({}, 0., 1.) { discr_mode_ = Discrete::Defined; }

            void discretize() override final {}

        private:
            double norm_() const { return domain_.max.value() - domain_.min.value() + 1.; }

        protected:
            ArgRequirements arg_requirements_() const override final
            {
                return {{&domain_.min, Codomain::integers()}, {&domain_.max, Codomain::integers()}};
            }

            double raw_function_(double) const override final { return 1. / norm_(); }
            double raw_primitive_(double x) const override final { return (x - domain_.min.value() + 1.) / norm_(); }
            double raw_log_func_(double) const override final { return -log(norm_()); }
            double raw_log_ratio_(double, double) const override final { return 0.; }
            double raw_quantile_(double p) const override final { return domain_.min.value() + std::ceil(p * norm_()) - 1.; }
            double raw_rand_(RndGen &prng) const override final
            {
                return std::uniform_int_distribution<long>(
                    static_cast<long>(domain_.min.value()), static_cast<long>(domain_.max.value()))(prng);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(UniformDiscrete_, Discrete)

        //////////////////////////////  Normal  //////////////////////////////

        class Normal_ final : public Cloneable<Normal_, Distribution>
        {
        public:
            static constexpr const char *NAME = "normal";
            Normal_() : Cloneable({{"mean", 0.}, {"sd", 1.}}) {}

            void discretize() override final
            {
                discr_mode_ = domain_.min.is_empty() ? Discrete::Center : Discrete::Left;
            }
            void truncate() override final
            {
                helper_init_standard_domain_();
                helper_init_distr_truncation_();
            }

        protected:
            ArgRequirements arg_requirements_() const override final { return positive_args_({1}); }

            double raw_function_(double x) const override final
            {
                const double sd = arg(1);
                const double z = (x - arg(0)) / sd;
                return M_1_SQRT_2PI * exp(-0.5 * z * z) / sd;
            }
            double raw_primitive_(double x) const override final
            {
                return 0.5 * std::erfc(M_SQRT1_2 * (arg(0) - x) / arg(1));
            }
            double raw_ccdf_(double x) const override final
            {
                return 0.5 * std::erfc(M_SQRT1_2 * (x - arg(0)) / arg(1));
            }
            double raw_log_func_(double x) const override final
            {
                const double sd = arg(1);
                const double z = (x - arg(0)) / sd;
                return -0.5 * z * z - M_LN_SQRT_2PI - log(sd);
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                const double sd = arg(1);
                return (0.5 * (xo * xo - xn * xn) + (xn - xo) * arg(0)) / (sd * sd);
            }
            double raw_quantile_(double p) const override final
            {
                return arg(0) + arg(1) * boost::math::erfc_inv(2. - 2. * p, no_promote_policy{}) * M_SQRT2;
            }
            double raw_rand_(RndGen &prng) const override final
            {
                return arg(0) + arg(1) * standard_normal(prng);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Normal_, Continuous)

        //////////////////////////////  Lognormal  //////////////////////////////

        /// Lognormal distribution.
        class Lognormal_ final : public Cloneable<Lognormal_, Distribution>
        {
        public:
            static constexpr const char *NAME = "lognormal";
            /// Build with parametrizations:
            /// * 0 -> meanlog-sdlog;
            /// * 1 -> mean-sd.
            Lognormal_() : Cloneable({{"meanlog", 0.}, {"sdlog", 1.}}, 0.)
            {
                alts_[0].name = "meanlog-sdlog";
                register_alt_("mean-sd", {{"mean", na::DOUBLE}, {"sd", na::DOUBLE}});
            }

            INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

            /// Finalize, then store the effective arguments of "mean-sd" when both arguments are constant.
            void finalize() override final
            {
                meanlog_ = sdlog_ = ArgCache{};
                Function::finalize();
                const bool is_stored = alt_idx_ == 1 && are_args_constant_({0, 1});
                meanlog_.store(is_stored, [this]
                               { return eff_params_().first; });
                sdlog_.store(is_stored, [this]
                             { return eff_params_().second; });
            }

        private:
            ArgCache meanlog_; /// Effective meanlog of "mean-sd", stored by finalize().
            ArgCache sdlog_;   /// Effective sdlog of "mean-sd", stored by finalize().

        protected:
            /// * meanlog-sdlog -> `sdlog` positive;
            /// * mean-sd -> both positive.
            ArgRequirements arg_requirements_() const override final
            {
                if (alt_idx_ == 0)
                    return positive_args_({1});
                return positive_args_({0, 1});
            }

            /// Effective (meanlog, sdlog) for the active parametrization.
            /// * "meanlog-sdlog" -> arguments;
            /// * "mean-sd" -> stored or computed.
            std::pair<double, double> eff_params_() const
            {
                switch (alt_idx_)
                {
                case 1: // "mean-sd": both read mean and sd, stored together
                {
                    if (sdlog_.is_stored())
                        return {meanlog_.value, sdlog_.value};
                    const double m = arg(0), r = arg(1) / m;
                    const double v = log1p(r * r);
                    return {log(m) - 0.5 * v, sqrt(v)};
                }
                default: // case 0: "meanlog-sdlog"
                    return {arg(0), arg(1)};
                }
            }

            double raw_function_(double x) const override final
            {
                if (x == 0.)
                    return 0.;

                const auto [meanlog, sdlog] = eff_params_();
                const double z = (log(x) - meanlog) / sdlog;
                return exp(-0.5 * z * z) * M_1_SQRT_2PI / (x * sdlog);
            }
            double raw_primitive_(double x) const override final
            {
                const auto [meanlog, sdlog] = eff_params_();
                return 0.5 * std::erfc(M_SQRT1_2 * (meanlog - log(x)) / sdlog);
            }
            double raw_ccdf_(double x) const override final
            {
                const auto [meanlog, sdlog] = eff_params_();
                return 0.5 * std::erfc(M_SQRT1_2 * (log(x) - meanlog) / sdlog);
            }
            double raw_log_func_(double x) const override final
            {
                if (x == 0.)
                    return -inf::DOUBLE;

                const auto [meanlog, sdlog] = eff_params_();
                const double z = (log(x) - meanlog) / sdlog;
                return -0.5 * z * z - log(x * sdlog) - M_LN_SQRT_2PI;
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                // the origin carries no density: the ratio saturates
                if (xo == 0. || xn == 0.)
                    return std::copysign(inf::DOUBLE, xn - xo);

                const auto [meanlog, sdlog] = eff_params_();
                const double lxo = log(xo);
                const double lxn = log(xn);
                return (0.5 * (lxo * lxo - lxn * lxn) + (lxn - lxo) * meanlog) / (sdlog * sdlog) + lxo - lxn;
            }
            double raw_quantile_(double p) const override final
            {
                const auto [meanlog, sdlog] = eff_params_();
                return exp(meanlog + sdlog * boost::math::erfc_inv(2. - 2. * p, no_promote_policy{}) * M_SQRT2);
            }
            double raw_rand_(RndGen &prng) const override final
            {
                const auto [meanlog, sdlog] = eff_params_();
                return exp(meanlog + sdlog * standard_normal(prng));
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Lognormal_, Continuous)

        //////////////////////////////  Beta  //////////////////////////////

        class Beta_ final : public Cloneable<Beta_, Distribution>
        {
        public:
            static constexpr const char *NAME = "beta";
            Beta_()
                : Cloneable({{"shape1", na::DOUBLE}, {"shape2", na::DOUBLE}}, 0., 1.)
            {
                domain_.below = Interval::Open;
                domain_.above = Interval::Open;
            }

            INIT_DISCR_REFUSED()
            void truncate() override final
            {
                const double sx = domain_.min.value();
                const double dx = domain_.max.value();
                if (sx < 0)
                    msg::stop(
                        "Invoked invalid min value " + str::format(sx) + " for " + str() + ".");
                if (dx > 1)
                    msg::stop(
                        "Invoked invalid max value " + str::format(dx) + " for " + str() + ".");
                if (sx > 0 || dx < 1)
                {
                    helper_init_standard_domain_(sx == 0, dx == 1);
                    helper_init_distr_truncation_();
                }
            }

        protected:
            ArgRequirements arg_requirements_() const override final { return positive_args_({0, 1}); }

            double raw_function_(double x) const override final
            {
                const double alpha = arg(0);
                const double beta = arg(1);
                return exp((alpha - 1.) * log(x) + (beta - 1.) * log1p(-x) - lbeta(alpha, beta));
            }
            double raw_primitive_(double x) const override final
            {
                return boost::math::ibeta(arg(0), arg(1), x, no_promote_policy{});
            }
            double raw_ccdf_(double x) const override final
            {
                return boost::math::ibetac(arg(0), arg(1), x, no_promote_policy{});
            }
            double raw_log_func_(double x) const override final
            {
                const double alpha = arg(0);
                const double beta = arg(1);
                return (alpha - 1.) * log(x) + (beta - 1.) * log1p(-x) - lbeta(alpha, beta);
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                return (arg(0) - 1.) * log(xn / xo) + (arg(1) - 1.) * log((1. - xn) / (1. - xo));
            }
            double raw_quantile_(double p) const override final
            {
                return boost::math::ibeta_inv(arg(0), arg(1), p, no_promote_policy{});
            }
            double raw_rand_(RndGen &prng) const override final
            {
                const double x = standard_gamma(prng, arg(0));
                return x / (x + standard_gamma(prng, arg(1)));
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Beta_, Continuous)

        //////////////////////////////  Gamma  //////////////////////////////

        class Gamma_ final : public Cloneable<Gamma_, Distribution>
        {
        public:
            static constexpr const char *NAME = "gamma";
            /// Build with parametrizations:
            /// * 0 -> shape-scale;
            /// * 1 -> mean-sd.
            Gamma_() : Cloneable({{"shape", na::DOUBLE}, {"scale", 1.}}, 0.)
            {
                alts_[0].name = "shape-scale";
                register_alt_("mean-sd", {{"mean", na::DOUBLE}, {"sd", na::DOUBLE}});
            }

            INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

            /// Finalize, then store each value when the arguments it reads are constant.
            /// * effective arguments of "mean-sd" -> both arguments;
            /// * log-gamma of the effective shape -> `shape`, or both arguments of "mean-sd";
            /// * log of the effective scale -> `scale`, or both arguments of "mean-sd".
            void finalize() override final
            {
                shape_ = scale_ = lgamma_shape_ = log_scale_ = ArgCache{};
                Function::finalize();
                const bool both = are_args_constant_({0, 1});
                shape_.store(alt_idx_ == 1 && both, [this]
                             { return eff_params_().first; });
                scale_.store(alt_idx_ == 1 && both, [this]
                             { return eff_params_().second; });
                lgamma_shape_.store(alt_idx_ == 1 ? both : are_args_constant_({0}), [this]
                                    { return lgamma(eff_params_().first); });
                log_scale_.store(alt_idx_ == 1 ? both : are_args_constant_({1}), [this]
                                 { return log(eff_params_().second); });
            }

        private:
            ArgCache shape_;        /// Effective shape of "mean-sd", stored by finalize().
            ArgCache scale_;        /// Effective scale of "mean-sd", stored by finalize().
            ArgCache lgamma_shape_; /// Log-gamma of the effective shape, stored by finalize().
            ArgCache log_scale_;    /// Log of the effective scale, stored by finalize().

        protected:
            /// Both arguments positive, in either parametrization.
            ArgRequirements arg_requirements_() const override final { return positive_args_({0, 1}); }

            /// Effective (shape, scale) for the active parametrization.
            /// * "shape-scale" -> arguments;
            /// * "mean-sd" -> stored or computed.
            std::pair<double, double> eff_params_() const
            {
                switch (alt_idx_)
                {
                case 1: // "mean-sd": both read mean and sd, stored together
                {
                    if (scale_.is_stored())
                        return {shape_.value, scale_.value};
                    const double sd = arg(1), r = arg(0) / sd;
                    return {r * r, sd / r};
                }
                default: // case 0: "shape-scale"
                    return {arg(0), arg(1)};
                }
            }

            double raw_function_(double x) const override final
            {
                if (x != 0.)
                    return exp(raw_log_func_(x));

                const auto [shape, scale] = eff_params_();
                return origin_pdf(shape, scale);
            }
            double raw_primitive_(double x) const override final
            {
                const auto [shape, scale] = eff_params_();
                return boost::math::gamma_p(shape, x / scale, no_promote_policy{});
            }
            double raw_ccdf_(double x) const override final
            {
                const auto [shape, scale] = eff_params_();
                return boost::math::gamma_q(shape, x / scale, no_promote_policy{});
            }
            double raw_log_func_(double x) const override final
            {
                const auto [shape, scale] = eff_params_();
                if (x == 0.)
                    return log(origin_pdf(shape, scale));
                return (shape - 1.) * log(x) - shape * log_scale_.read([s = scale]
                                                                       { return log(s); }) -
                       lgamma_shape_.read([a = shape]
                                          { return lgamma(a); }) -
                       x / scale;
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                // the origin density depends on the shape: difference of log-densities
                if (xo == 0. || xn == 0.)
                    return raw_log_func_(xn) - raw_log_func_(xo);

                const auto [shape, scale] = eff_params_();
                return (shape - 1.) * log(xn / xo) + (xo - xn) / scale;
            }
            double raw_quantile_(double p) const override final
            {
                const auto [shape, scale] = eff_params_();
                return scale * boost::math::gamma_p_inv(shape, p, no_promote_policy{});
            }
            double raw_rand_(RndGen &prng) const override final
            {
                const auto [shape, scale] = eff_params_();
                return scale * standard_gamma(prng, shape);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Gamma_, Continuous)

        //////////////////////////////  Rayleigh  //////////////////////////////

        class Rayleigh_ final : public Cloneable<Rayleigh_, Distribution>
        {
        public:
            static constexpr const char *NAME = "rayleigh";
            Rayleigh_() : Cloneable({{"scale", 1.}}, 0.) {}

            INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

        protected:
            ArgRequirements arg_requirements_() const override final { return positive_args_({0}); }

            double raw_function_(double x) const override final
            {
                const double scale = arg(0);
                const double scale2 = scale * scale;
                return x * exp(-0.5 * x * x / scale2) / scale2;
            }
            double raw_primitive_(double x) const override final
            {
                const double scale = arg(0);
                return -math::expm1(-0.5 * x * x / (scale * scale));
            }
            double raw_ccdf_(double x) const override final
            {
                const double scale = arg(0);
                return exp(-0.5 * x * x / (scale * scale));
            }
            double raw_log_func_(double x) const override final
            {
                const double scale = arg(0);
                const double scale2 = scale * scale;
                return log(x / scale2) - 0.5 * x * x / scale2;
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                const double scale = arg(0);
                return log(xn / xo) + 0.5 * (xo * xo - xn * xn) / (scale * scale);
            }
            double raw_quantile_(double p) const override final
            {
                return arg(0) * sqrt(-2. * log1p(-p));
            }
            double raw_rand_(RndGen &prng) const override final
            {
                return arg(0) * sqrt(-2 * log(standard_uniform(prng)));
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Rayleigh_, Continuous)

        //////////////////////////////  Weibull  //////////////////////////////

        class Weibull_ final : public Cloneable<Weibull_, Distribution>
        {
        public:
            static constexpr const char *NAME = "weibull";
            Weibull_() : Cloneable({{"shape", na::DOUBLE}, {"scale", 1.}}, 0.) {}

            INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

        protected:
            ArgRequirements arg_requirements_() const override final { return positive_args_({0, 1}); }

            double raw_function_(double x) const override final
            {
                const double shape = arg(0);
                const double scale = arg(1);
                if (x == 0.)
                    return origin_pdf(shape, scale);
                const double zk = pow(x / scale, shape);
                return shape / x * zk * exp(-zk);
            }
            double raw_primitive_(double x) const override final
            {
                return -math::expm1(-pow(x / arg(1), arg(0)));
            }
            double raw_ccdf_(double x) const override final
            {
                return exp(-pow(x / arg(1), arg(0)));
            }
            double raw_log_func_(double x) const override final
            {
                const double shape = arg(0);
                const double scale = arg(1);
                if (x == 0.)
                    return log(origin_pdf(shape, scale));
                const double z = x / scale;
                return log(shape / scale) + (shape - 1.) * log(z) - pow(z, shape);
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                // the origin density depends on the shape: difference of log-densities
                if (xo == 0. || xn == 0.)
                    return raw_log_func_(xn) - raw_log_func_(xo);

                const double shape = arg(0);
                const double scale = arg(1);
                return (shape - 1.) * log(xn / xo) + pow(xo / scale, shape) - pow(xn / scale, shape);
            }
            double raw_quantile_(double p) const override final
            {
                return arg(1) * pow(-log1p(-p), 1. / arg(0));
            }
            double raw_rand_(RndGen &prng) const override final
            {
                return arg(1) * pow(-log(standard_uniform(prng)), 1.0 / arg(0));
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Weibull_, Continuous)

        //////////////////////////////  Exponential  //////////////////////////////

        class Exponential_ final : public Cloneable<Exponential_, Distribution>
        {
        public:
            static constexpr const char *NAME = "exponential";
            Exponential_() : Cloneable({{"rate", 1.}}, 0.) {}

            INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

        protected:
            ArgRequirements arg_requirements_() const override final { return positive_args_({0}); }

            double raw_function_(double x) const override final
            {
                const double rate = arg(0);
                return rate * exp(-rate * x);
            }
            double raw_primitive_(double x) const override final
            {
                return -math::expm1(-arg(0) * x);
            }
            double raw_ccdf_(double x) const override final
            {
                return exp(-arg(0) * x);
            }
            double raw_log_func_(double x) const override final
            {
                const double rate = arg(0);
                return log(rate) - rate * x;
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                return arg(0) * (xo - xn);
            }
            double raw_quantile_(double p) const override final
            {
                return -log1p(-p) / arg(0);
            }
            double raw_rand_(RndGen &prng) const override final
            {
                return -log(standard_uniform(prng)) / arg(0);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Exponential_, Continuous)

        //////////////////////////////  Poisson  //////////////////////////////

        class Poisson_ final : public Cloneable<Poisson_, Distribution>
        {
        public:
            static constexpr const char *NAME = "poisson";
            Poisson_() : Cloneable({{"lambda", na::DOUBLE}}, 0.)
            {
                discr_mode_ = Discrete::Defined;
            }

            INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

        protected:
            ArgRequirements arg_requirements_() const override final { return positive_args_({0}); }

            double raw_function_(double x) const override final
            {
                return exp(raw_log_func_(x));
            }
            double raw_primitive_(double x) const override final
            {
                return boost::math::gamma_q(x + 1., arg(0), no_promote_policy{});
            }
            double raw_log_func_(double x) const override final
            {
                const double lambda = arg(0);
                return x * log(lambda) - lambda - lgamma(x + 1.);
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                return (xn - xo) * log(arg(0)) + lgamma(xo + 1.) - lgamma(xn + 1.);
            }
            double raw_quantile_(double p) const override final
            {
                return boost::math::quantile(
                    boost::math::poisson_distribution<double, discrete_up_policy>(arg(0)), p);
            }
            double raw_rand_(RndGen &prng) const override final
            {
                return std::poisson_distribution<unsigned>(arg(0))(prng);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Poisson_, Discrete)

        //////////////////////////////  Binomial  //////////////////////////////

        class Binomial_ final : public Cloneable<Binomial_, Distribution>
        {
        public:
            static constexpr const char *NAME = "binomial";
            Binomial_()
                : Cloneable({{"size", na::DOUBLE}, {"prob", na::DOUBLE}}, 0.)
            {
                discr_mode_ = Discrete::Defined;
                domain_.max.init_pointer(arguments_(0).pointer());
                domain_.below = Interval::Closed;
                domain_.above = Interval::Closed;
            }

            void discretize() override final {}
            void truncate() override final
            {
                const double sx = domain_.min.value();
                if (sx < 0.)
                    msg::stop(
                        "Invoked invalid min value " + str::format(sx) +
                        " for " + str() + ".");
                if (sx > 0. || arguments_(0).pointer() != max().pointer())
                {
                    const double dx = domain_.max.value();
                    if (dx > arg(0))
                        msg::stop(
                            "Invoked invalid max value " + str::format(dx) +
                            " for " + str() + " with size " + str::format(arg(0)) + ".");
                    helper_init_standard_domain_();
                    helper_init_distr_truncation_();
                }
            }

        protected:
            /// `size` counts trials.
            ArgRequirements arg_requirements_() const override final
            {
                return {{&arguments_(0), Codomain::non_negative(Lattice::Integers)},
                        {&arguments_(1), Codomain::unit(Interval::Open, Interval::Open)}};
            }

            double raw_function_(double x) const override final
            {
                const double size = arg(0);
                const double prob = arg(1);
                return exp(log_binom_coeff(size, x) + x * log(prob) + (size - x) * log1p(-prob));
            }
            double raw_primitive_(double x) const override final
            {
                const double size = arg(0);
                if (x > size)
                    return 1.;
                return boost::math::ibeta(size - x, x + 1, 1. - arg(1), no_promote_policy{});
            }
            double raw_log_func_(double x) const override final
            {
                const double size = arg(0);
                const double prob = arg(1);
                return log_binom_coeff(size, x) + x * log(prob) + (size - x) * log1p(-prob);
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                const double prob = arg(1);
                return logratio_binom_coeff(arg(0), xo, xn) + (xn - xo) * log(prob / (1. - prob));
            }
            double raw_quantile_(double p) const override final
            {
                return boost::math::quantile(
                    boost::math::binomial_distribution<double, discrete_up_policy>(arg(0), arg(1)), p);
            }
            double raw_rand_(RndGen &prng) const override final
            {
                return std::binomial_distribution<unsigned>(arg(0), arg(1))(prng);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Binomial_, Discrete)

        //////////////////////////////  NegativeBinomial  //////////////////////////////

        class NegativeBinomial_ final : public Cloneable<NegativeBinomial_, Distribution>
        {
        public:
            static constexpr const char *NAME = "negative_binomial";
            NegativeBinomial_()
                : Cloneable({{"size", na::DOUBLE}, {"prob", na::DOUBLE}}, 0.)
            {
                discr_mode_ = Discrete::Defined;
            }

            INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

        protected:
            /// `size` counts successes: the binomial coefficient truncates it to an integer.
            ArgRequirements arg_requirements_() const override final
            {
                return {{&arguments_(0), Codomain::positive(Lattice::Integers)},
                        {&arguments_(1), Codomain::unit(Interval::Open, Interval::Open)}};
            }

            double raw_function_(double x) const override final
            {
                const double size = arg(0);
                const double prob = arg(1);
                return exp(log_binom_coeff(x + size - 1., x) + x * log1p(-prob) + size * log(prob));
            }
            double raw_primitive_(double x) const override final
            {
                return boost::math::ibeta(arg(0), x + 1., arg(1), no_promote_policy{});
            }
            double raw_log_func_(double x) const override final
            {
                const double size = arg(0);
                const double prob = arg(1);
                return log_binom_coeff(x + size - 1., x) + x * log1p(-prob) + size * log(prob);
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                const double size = arg(0);
                const double prob = arg(1);
                return lgamma(xn + size) + lgamma(xo + 1) - lgamma(xo + size) - lgamma(xn + 1) + (xn - xo) * log1p(-prob);
            }
            double raw_quantile_(double p) const override final
            {
                return boost::math::quantile(
                    boost::math::negative_binomial_distribution<double, discrete_up_policy>(arg(0), arg(1)), p);
            }
            double raw_rand_(RndGen &prng) const override final
            {
                return std::negative_binomial_distribution<unsigned>(arg(0), arg(1))(prng);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(NegativeBinomial_, Discrete)

        //////////////////////////////  Geometric  //////////////////////////////

        class Geometric_ final : public Cloneable<Geometric_, Distribution>
        {
        public:
            static constexpr const char *NAME = "geometric";
            Geometric_() : Cloneable({{"prob", na::DOUBLE}}, 0.)
            {
                discr_mode_ = Discrete::Defined;
            }

            INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

        private:
            /// Shrink of `p` keeping an accurate cumulative atom on its own quantile.
            static constexpr double QUANTILE_SHRINK = 1. - 2. * std::numeric_limits<double>::epsilon();

        protected:
            ArgRequirements arg_requirements_() const override final
            {
                return {{&arguments_(0), Codomain::unit(Interval::Open, Interval::Open)}};
            }

            double raw_function_(double x) const override final
            {
                const double prob = arg(0);
                return prob * exp(x * log1p(-prob));
            }
            double raw_primitive_(double x) const override final
            {
                return -math::expm1((x + 1.) * log1p(-arg(0)));
            }
            double raw_log_func_(double x) const override final
            {
                const double prob = arg(0);
                return log(prob) + x * log1p(-prob);
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                return (xn - xo) * log1p(-arg(0));
            }
            double raw_quantile_(double p) const override final
            {
                // closed form overshoots by one where the cumulative rounds below `p`: step back
                const double k = std::ceil(log1p(-p * QUANTILE_SHRINK) / log1p(-arg(0)) - 1.);
                return k > 0. && raw_primitive_(k - 1.) >= p ? k - 1. : k;
            }
            double raw_rand_(RndGen &prng) const override final
            {
                return std::geometric_distribution<unsigned>(arg(0))(prng);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Geometric_, Discrete)

        //////////////////////////////  Categorical  //////////////////////////////

        class Categorical_ final : public Cloneable<Categorical_, Distribution>
        {
        public:
            static constexpr const char *NAME = "categorical";
            Categorical_()
            {
                discr_mode_ = Discrete::Defined;
            }

            void discretize() override final {}
            void truncate() override final
            {
                const size_t cases = n_args();
                if (cases == 0)
                    msg::stop(
                        "Invoked " + str() + " with no weights.\nProvide at least one weight.");

                if (norm_() <= 0.)
                    msg::stop(
                        "Invoked " + str() + " with weights summing to non-strictly positive value.");

                domain_.min.init_numeric(0.);
                domain_.max.init_numeric(static_cast<double>(cases - 1));
                domain_.below = Interval::Closed;
                domain_.above = Interval::Closed;
            }

        private:
            double norm_() const
            {
                double z = 0.;
                for (size_t i = 0; i < n_args(); ++i)
                    z += arg(i);
                return z;
            }

        protected:
            /// Every weight non-negative.
            ArgRequirements arg_requirements_() const override final
            {
                ArgRequirements reqs;
                for (const Attribute &attr : arguments_())
                    reqs.emplace_back(&attr, Codomain::non_negative());
                return reqs;
            }

            Attribute &find_arg_(std::string_view key) override final { return append_arg_(key); }

            double raw_function_(double x) const override final
            {
                return arg(static_cast<size_t>(x)) / norm_();
            }
            double raw_log_func_(double x) const override final
            {
                return log(arg(static_cast<size_t>(x)) / norm_());
            }
            double raw_primitive_(double x) const override final
            {
                const size_t k = static_cast<size_t>(x);
                double acc = 0.;
                for (size_t i = 0; i <= k; ++i)
                    acc += arg(i);
                return acc / norm_();
            }
            double raw_log_ratio_(double xo, double xn) const override final
            {
                return log(arg(static_cast<size_t>(xn)) / arg(static_cast<size_t>(xo)));
            }
            double raw_quantile_(double p) const override final
            {
                const size_t cases = n_args();
                const double target = p * norm_();
                size_t i = 0;
                double acc = arg(0);
                while (i < cases && acc < target)
                    acc += arg(++i);
                return static_cast<double>(i);
            }
            double raw_rand_(RndGen &prng) const override final
            {
                const size_t cases = n_args();
                const double target = standard_uniform(prng) * norm_();
                size_t i = 0;
                double acc = arg(0);
                while (i < cases && acc < target)
                    acc += arg(++i);
                return static_cast<double>(i);
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Categorical_, Discrete)

        //////////////////////////////  Ferretti20  //////////////////////////////

        /// Infectivity profile as the distribution of the time from symptom onset to transmission.
        /// Based on @cite ferretti_timing_2020
        ///
        /// * `t_ost`: time of infection of the infectee minus time of onset of the infector;
        /// * `t_ost` negative -> pre-symptomatic transmission;
        /// * skew-logistic distribution, generalized logistic type I;
        /// * pre-symptomatic values rescaled according to the incubation periods.
        ///
        /// Arguments, in order:
        /// * `shape`, `scale`;
        /// * `c_loc`: custom location;
        /// * `mean_incub`: mean incubation period;
        /// * `incub_period`: individual incubation period.
        class Ferretti20_ final : public Cloneable<Ferretti20_, Distribution>
        {
        public:
            static constexpr const char *NAME = "ferretti20";
            Ferretti20_()
                : Cloneable({{"shape", 5.85}, {"scale", 1.85}, {"c_loc", -4.0}, {"mean_incub", 5.42}, {"incub_period", na::DOUBLE}}, 0.) {}

            INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

            /// Finalize, then store each onset term when the arguments it reads are constant.
            void finalize() override final
            {
                pema_ = pemTSa_ = ArgCache{};
                Function::finalize();
                pema_.store(are_args_constant_({0, 1, 2}), [this]
                            { return onset_term_(arg(0), arg(1), arg(2)); });
                pemTSa_.store(are_args_constant_({0, 1, 2, 3}), [this]
                              { return onset_term_(arg(0), arg(1), arg(2) + arg(3)); });
            }

        private:
            ArgCache pema_;   /// Onset term at `c_loc`, stored by finalize().
            ArgCache pemTSa_; /// Onset term at `c_loc + mean_incub`, stored by finalize().

            /// Onset term (1 + exp(`m` / `s`))^-`a`.
            static double onset_term_(double a, double s, double m) { return pow(1. + exp(m / s), -a); }
            /// Onset terms at `c_loc` = `m` and `c_loc + mean_incub` = `m + T`, stored or computed.
            std::pair<double, double> onset_terms_(double a, double s, double m, double T) const
            {
                return {pema_.read([=]
                                   { return onset_term_(a, s, m); }),
                        pemTSa_.read([=]
                                     { return onset_term_(a, s, m + T); })};
            }

        protected:
            /// Every argument positive but `c_loc`.
            ArgRequirements arg_requirements_() const override final
            {
                return positive_args_({0, 1, 3, 4});
            }

            double raw_rand_(RndGen & /*prng*/) const override final { raise_undefined_call_err_("rand"); }
            double raw_log_ratio_(double /*xo*/, double /*xn*/) const override final { raise_undefined_call_err_("log_ratio"); }
            double raw_quantile_(double /*p*/) const override final { raise_undefined_call_err_("quantile"); }
            double raw_log_func_(double /*x*/) const override final { raise_undefined_call_err_("log_pdf"); }
            double raw_function_(double x) const override final
            {
                const double a = arg(0), s = arg(1), m = arg(2), T = arg(3), i = arg(4);
                const double tOST = x - i; // time from onset to transmission
                const double iT = i / T;   // re-scaled incubation period
                const auto [pema, pemTSa] = onset_terms_(a, s, m, T);

                const double norm_pdf = (s / a) * (1. + (iT - 1.) * pema - iT * pemTSa);
                const double emOSTs = exp((m - (tOST < 0. ? tOST / iT : tOST)) / s);
                return emOSTs * pow(1. + emOSTs, -a - 1.) / norm_pdf;
            }
            double raw_primitive_(double x) const override final
            {
                const double a = arg(0), s = arg(1), m = arg(2), T = arg(3), i = arg(4);
                const double tOST = x - i; // time from onset to transmission
                const double iT = i / T;   // re-scaled incubation period
                const auto [pema, pemTSa] = onset_terms_(a, s, m, T);
                const double norm_cdf = 1. + (iT - 1.) * pema - iT * pemTSa;

                if (tOST < 0) // integrate from -incub_period to tOST
                    return iT * (pow(1. + exp((m - tOST / iT) / s), -a) - pemTSa) / norm_cdf;

                // integrate from -incub_period to 0
                const double F0 = iT * (pema - pemTSa);
                // integrate from 0 to tOST
                return (F0 + pow(1. + exp((m - tOST) / s), -a) - pema) / norm_cdf;
            }
            double raw_ccdf_(double x) const override final
            {
                const double a = arg(0), s = arg(1), m = arg(2), T = arg(3), i = arg(4);
                const double tOST = x - i; // time from onset to transmission
                const double iT = i / T;   // re-scaled incubation period
                const auto [pema, pemTSa] = onset_terms_(a, s, m, T);
                const double norm_cdf = 1. + (iT - 1.) * pema - iT * pemTSa;

                if (tOST < 0) // integrate from tOST to infinity
                    return (1. + (iT - 1.) * pema - iT * pow(1. + exp((m - tOST / iT) / s), -a)) / norm_cdf;

                // integrate from tOST to infinity, 1 - (1 + exp((m - tOST) / s))^-a without cancellation
                return -math::expm1(-a * log1p(exp((m - tOST) / s))) / norm_cdf;
            }
        };
        HOUSEPI_REGISTER_FUNCTION(Ferretti20_, InfecProfile)
    }
}
