/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/random/beta_distribution.hpp>

/// Block boost long double policy.
using no_promote_policy = boost::math::policies::policy<
    boost::math::policies::promote_float<false>,
    boost::math::policies::promote_double<false>>;

/// Discrete quantile policy matching R's convention: min{k : P(X<=k) >= p}.
/// Boost's default 'integer_round_outwards' rounds toward the distribution tails
/// which differs from R for p < 0.5. 'integer_round_up' matches R exactly.
using discrete_up_policy = boost::math::policies::policy<
    boost::math::policies::promote_float<false>,
    boost::math::policies::promote_double<false>,
    boost::math::policies::discrete_quantile<boost::math::policies::integer_round_up>>;

#include "tools/messages.h"
#include "tools/string.h"

#include "distribution.h"
#include "registry.h"

namespace housepi
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////// *** DistributionRegistry : class ***  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    std::unique_ptr<Function> DistributionRegistry::create(std::string_view opt, bool is_prim) const
    {
        auto it = factory_map_.find(opt);
        if (it == factory_map_.end())
            return nullptr;

        auto output = it->second();
        if (is_prim)
            output->set_eval_as_cdf();
        else
            output->set_eval_as_pdf();

        return output;
    }

    std::unique_ptr<Distribution> compose_distr_registries(std::string_view opt,
                                                           const std::string &err_type,
                                                           const std::vector<std::string> &registries,
                                                           bool is_cumul)
    {
        if (registries.empty())
            msg::stop(
                "No registries were provided for " + msg::val(err_type) + ".");

        // FIXME: this should be an error
        if (na::is(opt))
            return nullptr;

        std::vector<const DistributionRegistry *> registry_ptrs = {};
        for (const std::string &r : registries)
            if (r == "continuos")
                registry_ptrs.push_back(&ContinuosDistribRegistry::instance());
            else if (r == "discrete")
                registry_ptrs.push_back(&DiscreteDistribRegistry::instance());
            else if (r == "infecprofile")
                registry_ptrs.push_back(&InfecProfileRegistry::instance());
            else
                msg::stop(
                    "Unknown " + err_type + " registry " + msg::val(r) + ".");

        std::vector<std::string> options = {};
        for (const DistributionRegistry *r : registry_ptrs)
        {
            if (r->count(opt))
                return std::unique_ptr<Distribution>(
                    static_cast<Distribution *>(r->create(opt, is_cumul).release()));
            const std::vector<std::string> &y = r->get_options();
            options.reserve(options.size() + y.size());
            options.insert(options.end(), y.begin(), y.end());
        }

        msg::stop(
            "Unknown " + err_type + " " + msg::val(std::string(opt)) + "." +
            str::find_closest(opt, options));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  *** Distribution : implementations ***  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    using namespace housepi;
    using std::exp;
    using std::lgamma;
    using std::log;
    using std::pow;
    using std::sqrt;

#ifndef M_1_SQRT_2PI
    constexpr double M_1_SQRT_2PI = 0.398942280401432677939946059934; /// 1/sqrt(2pi)
#endif

#ifndef M_SQRT2
    constexpr double M_SQRT2 = 1.41421356237309504880168872420969808; /// sqrt(2)
#endif

#ifndef M_SQRT1_2
    constexpr double M_SQRT1_2 = 0.707106781186547524400844362105; /// 1/sqrt(2)
#endif

#ifndef M_LN_SQRT_2PI
    constexpr double M_LN_SQRT_2PI = 0.918938533204672741780329736406; /// log(sqrt(2pi))
#endif

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

// Macro for left discretization when min = 0 (close).
#define INIT_DISCR_LEFT_MIN_ZERO_(is_open)       \
    void init_discretization() override final    \
    {                                            \
        discr_mode_ = Mode::Discrete::Left;      \
    }                                            \
    void init_truncation() override final        \
    {                                            \
        if (helper_init_range_min_(0., is_open)) \
            helper_init_distr_truncation_();     \
    }
#define INIT_DISCR_LEFT_MIN_ZERO_CLOSE() \
    INIT_DISCR_LEFT_MIN_ZERO_(false)
#define INIT_DISCR_LEFT_MIN_ZERO_OPEN() \
    INIT_DISCR_LEFT_MIN_ZERO_(true)

// Macro for discrete defined when min = 0 (close).
#define INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()      \
    void init_discretization() override final {} \
    void init_truncation() override final        \
    {                                            \
        if (helper_init_range_min_(0., false))   \
            helper_init_distr_truncation_();     \
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  *** Uniform Distribution : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Uniform_ final : public Distribution
    {
    public:
        Uniform_()
            : Distribution({}, 0., 1.),
              unif_gen(0., 1.) {}
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Uniform_)
        HOUSEPI_OVERRIDE_NAME_METHOD(uniform)

        /// UPDATE: should understand when discrete defined...
        void init_discretization() override final { discr_mode_ = Mode::Discrete::Center; }

    protected:
        mutable std::uniform_real_distribution<double> unif_gen;

        double raw_pdf_(double) const override final { return 1. / (upper_.value() - lower_.value()); }
        double raw_cdf_(double x) const override final { return (x - lower_.value()) / (upper_.value() - lower_.value()); }
        double raw_log_pdf_(double) const override final { return -std::log(upper_.value() - lower_.value()); }
        double raw_log_ratio_(double, double) const override final { return 0.; }
        double raw_rand_(rnd_gen &prng) const override final
        {
            const double r = unif_gen(prng);
            return lower_.value() * (1. - r) + upper_.value() * r;
        }
        double raw_quantile_(double p) const override final { return lower_.value() * (1. - p) + upper_.value() * p; }
    };
    HOUSEPI_REGISTER_FUNCTION(Uniform_, ContinuosDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Normal Distribution : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Normal_ final : public Distribution
    {
    public:
        Normal_()
            : Distribution({{"mean", 0.}, {"sd", 1.}}),
              norm_gen(0., 1.) {}
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Normal_)
        HOUSEPI_OVERRIDE_NAME_METHOD(normal)

        void init_discretization() override final
        {
            discr_mode_ = lower_.is_empty() ? Mode::Discrete::Center
                                            : Mode::Discrete::Left;
        }
        void init_truncation() override final
        {
            helper_init_standard_range_();
            helper_init_distr_truncation_();
        }

    protected:
        mutable std::normal_distribution<double> norm_gen;

        double raw_pdf_(double x) const override final
        {
            const double sd = arg(1);
            const double z = (x - arg(0)) / sd;
            return M_1_SQRT_2PI * exp(-0.5 * z * z) / sd;
        }
        double raw_cdf_(double x) const override final
        {
            return 0.5 * std::erfc(M_SQRT1_2 * (arg(0) - x) / arg(1));
        }
        double raw_ccdf_(double x) const override final
        {
            return 0.5 * std::erfc(M_SQRT1_2 * (x - arg(0)) / arg(1));
        }
        double raw_log_pdf_(double x) const override final
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
        double raw_rand_(rnd_gen &prng) const override final
        {
            return arg(0) + arg(1) * norm_gen(prng);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Normal_, ContinuosDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  *** Lognormal Distribution : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Lognormal Distribution derived class.
    class Lognormal_ final : public Distribution
    {
    public:
        Lognormal_()
            : Distribution({{"meanlog", 0.}, {"sdlog", 1.}}, 0.),
              norm_gen(0., 1.) { below_mode_ = Mode::Interval::Open; }
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Lognormal_)
        HOUSEPI_OVERRIDE_NAME_METHOD(lognormal)

        INIT_DISCR_LEFT_MIN_ZERO_OPEN()

    protected:
        mutable std::normal_distribution<double> norm_gen;

        double raw_pdf_(double x) const override final
        {
            const double sdlog = arg(1);
            const double z = (log(x) - arg(0)) / sdlog;
            return exp(-0.5 * z * z) * M_1_SQRT_2PI / (x * sdlog);
        }
        double raw_cdf_(double x) const override final
        {
            return 0.5 * std::erfc(M_SQRT1_2 * (arg(0) - log(x)) / arg(1));
        }
        double raw_ccdf_(double x) const override final
        {
            return 0.5 * std::erfc(M_SQRT1_2 * (log(x) - arg(0)) / arg(1));
        }
        double raw_log_pdf_(double x) const override final
        {
            const double sdlog = arg(1);
            const double z = (log(x) - arg(0)) / sdlog;
            return -0.5 * z * z - log(x * sdlog) - M_LN_SQRT_2PI;
        }
        double raw_log_ratio_(double xo, double xn) const override final
        {
            const double lxo = log(xo);
            const double lxn = log(xn);
            const double sdlog = arg(1);
            return (0.5 * (lxo * lxo - lxn * lxn) + (lxn - lxo) * arg(0)) / (sdlog * sdlog) + lxo - lxn;
        }
        double raw_quantile_(double p) const override final
        {
            return exp(arg(0) + arg(1) * boost::math::erfc_inv(2. - 2. * p, no_promote_policy{}) * M_SQRT2);
        }
        double raw_rand_(rnd_gen &prng) const override final
        {
            return exp(arg(0) + arg(1) * norm_gen(prng));
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Lognormal_, ContinuosDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** Beta Distribution : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Beta_ final : public Distribution
    {
    public:
        Beta_()
            : Distribution({{"shape1", na::DOUBLE}, {"shape2", na::DOUBLE}},
                           0., 1.)
        {
            below_mode_ = Mode::Interval::Open;
            above_mode_ = Mode::Interval::Open;
        }
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Beta_)
        HOUSEPI_OVERRIDE_NAME_METHOD(beta)

        void init_discretization() override final
        {
            msg::stop(
                "Invoked invalid option discrete for " + err_name() + ".");
        }
        void init_truncation() override final
        {
            const double sx = lower_.value();
            const double dx = upper_.value();
            if (sx < 0)
                msg::stop(
                    "Invoked invalid min value " + str::format(sx) + " for " + err_name() + ".");
            if (dx > 1)
                msg::stop(
                    "Invoked invalid max value " + str::format(dx) + " for " + err_name() + ".");
            if (sx > 0 || dx < 1)
            {
                helper_init_standard_range_(sx == 0, dx == 1);
                helper_init_distr_truncation_();
            }
        }

    protected:
        double raw_pdf_(double x) const override final
        {
            const double alpha = arg(0);
            const double beta = arg(1);
            return exp((alpha - 1.) * log(x) + (beta - 1.) * log1p(-x) - lbeta(alpha, beta));
        }
        double raw_cdf_(double x) const override final
        {
            return boost::math::ibeta(arg(0), arg(1), x, no_promote_policy{});
        }
        double raw_log_pdf_(double x) const override final
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
        double raw_rand_(rnd_gen &prng) const override final
        {
            return boost::random::beta_distribution<double>(arg(0), arg(1))(prng);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Beta_, ContinuosDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** Gamma Distribution : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Gamma_ final : public Distribution
    {
    public:
        /// Alternative arametrizations:
        /// 0. "shape-scale"
        /// 1. "mean-var"
        Gamma_()
            : Distribution({{"shape", na::DOUBLE}, {"scale", 1.}}, 0.)
        {
            alts_[0].name = "shape-scale";
            register_alt_("mean-var", {{"mean", na::DOUBLE}, {"var", na::DOUBLE}});
        }
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Gamma_)
        HOUSEPI_OVERRIDE_NAME_METHOD(gamma)

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        /// Effective (shape, scale) for the active parametrization.
        std::pair<double, double> eff_params_() const
        {
            switch (alt_idx_)
            {
            case 1: // "mean-var"
            {
                const double m = arg(0), v = arg(1);
                return {m * m / v, v / m};
            }
            default: // case 0: "shape-scale"
                return {arg(0), arg(1)};
            }
        }

        double raw_pdf_(double x) const override final
        {
            if (x != 0.)
                return exp(raw_log_pdf_(x));

            const auto [shape, scale] = eff_params_();
            if (shape < 1.)
                return std::numeric_limits<double>::infinity();
            else if (shape == 1.)
                return exp(-log(scale) - lgamma(1));
            else
                return 0.;
        }
        double raw_cdf_(double x) const override final
        {
            const auto [shape, scale] = eff_params_();
            return boost::math::gamma_p(shape, x / scale, no_promote_policy{});
        }
        double raw_ccdf_(double x) const override final
        {
            const auto [shape, scale] = eff_params_();
            return boost::math::gamma_q(shape, x / scale, no_promote_policy{});
        }
        double raw_log_pdf_(double x) const override final
        {
            const auto [shape, scale] = eff_params_();
            if (x == 0.)
            {
                if (shape < 1.)
                    return std::numeric_limits<double>::infinity();
                if (shape == 1.)
                    return -std::log(scale);
                return -std::numeric_limits<double>::infinity();
            }
            const double t = x / scale;
            return shape * log(t) - log(x) - lgamma(shape) - t;
        }
        double raw_log_ratio_(double xo, double xn) const override final
        {
            const auto [shape, scale] = eff_params_();
            return (shape - 1.) * log(xn / xo) + (xo - xn) / scale;
        }
        double raw_quantile_(double p) const override final
        {
            const auto [shape, scale] = eff_params_();
            return scale * boost::math::gamma_p_inv(shape, p, no_promote_policy{});
        }
        double raw_rand_(rnd_gen &prng) const override final
        {
            const auto [shape, scale] = eff_params_();
            return std::gamma_distribution<double>(shape, scale)(prng);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Gamma_, ContinuosDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Rayleigh Distribution : class ***  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Rayleigh_ final : public Distribution
    {
    public:
        Rayleigh_()
            : Distribution({{"scale", 1.}}, 0.),
              unif_gen(0., 1.) {}
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Rayleigh_)
        HOUSEPI_OVERRIDE_NAME_METHOD(rayleigh)

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        mutable std::uniform_real_distribution<double> unif_gen;

        double raw_pdf_(double x) const override final
        {
            const double scale = arg(0);
            const double scale2 = scale * scale;
            return x * exp(-0.5 * x * x / scale2) / scale2;
        }
        double raw_cdf_(double x) const override final
        {
            return 1. - raw_ccdf_(x);
        }
        double raw_ccdf_(double x) const override final
        {
            const double scale = arg(0);
            return exp(-0.5 * x * x / (scale * scale));
        }
        double raw_log_pdf_(double x) const override final
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
        double raw_rand_(rnd_gen &prng) const override final
        {
            return arg(0) * sqrt(-2 * log(unif_gen(prng)));
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Rayleigh_, ContinuosDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Weibull Distribution : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Weibull_ final : public Distribution
    {
    public:
        Weibull_()
            : Distribution({{"shape", na::DOUBLE}, {"scale", 1.}}, 0.),
              unif_gen(0., 1.) {}
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Weibull_)
        HOUSEPI_OVERRIDE_NAME_METHOD(weibull)

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        mutable std::uniform_real_distribution<double> unif_gen;

        double raw_pdf_(double x) const override final
        {
            const double shape = arg(0);
            const double scale = arg(1);
            const double z = x / scale;
            return (shape / scale) * pow(z, shape - 1.) * exp(-pow(z, shape));
        }
        double raw_cdf_(double x) const override final
        {
            return 1. - raw_ccdf_(x);
        }
        double raw_ccdf_(double x) const override final
        {
            return exp(-pow(x / arg(1), arg(0)));
        }
        double raw_log_pdf_(double x) const override final
        {
            const double shape = arg(0);
            const double scale = arg(1);
            const double z = x / scale;
            return log(shape / scale) + (shape - 1.) * log(z) - pow(z, shape);
        }
        double raw_log_ratio_(double xo, double xn) const override final
        {
            const double shape = arg(0);
            const double scale = arg(1);
            return (shape - 1.) * log(xn / xo) + pow(xo / scale, shape) - pow(xn / scale, shape);
        }
        double raw_quantile_(double p) const override final
        {
            return arg(1) * pow(-log1p(-p), 1. / arg(0));
        }
        double raw_rand_(rnd_gen &prng) const override final
        {
            return arg(1) * pow(-log(unif_gen(prng)), 1.0 / arg(0));
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Weibull_, ContinuosDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  *** Exponential Distribution : class ***  //////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Exponential_ final : public Distribution
    {
    public:
        Exponential_()
            : Distribution({{"rate", 1.}}, 0.),
              unif_gen(0., 1.) {}
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Exponential_)
        HOUSEPI_OVERRIDE_NAME_METHOD(exponential)

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        mutable std::uniform_real_distribution<double> unif_gen;

        double raw_pdf_(double x) const override final
        {
            const double rate = arg(0);
            return rate * exp(-rate * x);
        }
        double raw_cdf_(double x) const override final
        {
            return 1. - raw_ccdf_(x);
        }
        double raw_ccdf_(double x) const override final
        {
            return exp(-arg(0) * x);
        }
        double raw_log_pdf_(double x) const override final
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
        double raw_rand_(rnd_gen &prng) const override final
        {
            return -log(unif_gen(prng)) / arg(0);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Exponential_, ContinuosDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** Poisson Distribution : class ***  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Poisson_ final : public Distribution
    {
    public:
        Poisson_()
            : Distribution({{"lambda", na::DOUBLE}}, 0.) { discr_mode_ = Mode::Discrete::Defined; }
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Poisson_)
        HOUSEPI_OVERRIDE_NAME_METHOD(poisson)

        INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            return exp(raw_log_pdf_(x));
        }
        double raw_cdf_(double x) const override final
        {
            return 1. - boost::math::gamma_p(x + 1., arg(0), no_promote_policy{});
        }
        double raw_log_pdf_(double x) const override final
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
        double raw_rand_(rnd_gen &prng) const override final
        {
            return std::poisson_distribution<unsigned>(arg(0))(prng);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Poisson_, DiscreteDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** Binomial Distribution : class ***  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Binomial_ final : public Distribution
    {
    public:
        Binomial_()
            : Distribution({{"size", na::DOUBLE}, {"prob", na::DOUBLE}}, 0.)
        {
            discr_mode_ = Mode::Discrete::Defined;
            upper_.init_pointer(arguments_(0).pointer());
            below_mode_ = Mode::Interval::Closed;
            above_mode_ = Mode::Interval::Closed;
        }
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Binomial_)
        HOUSEPI_OVERRIDE_NAME_METHOD(binomial)

        void init_discretization() override final {}
        void init_truncation() override final
        {
            const double sx = lower_.value();
            if (sx < 0.)
                msg::stop(
                    "Invoked invalid min value " + str::format(sx) +
                    " for " + err_name() + ".");
            if (sx > 0. || arguments_(0).pointer() != upper().pointer())
            {
                const double dx = lower_.value();
                if (dx > arg(0))
                    msg::stop(
                        "Invoked invalid max value " + str::format(dx) +
                        " for " + err_name() + " with size " + std::to_string(arg(0)) + ".");
                helper_init_standard_range_();
                helper_init_distr_truncation_();
            }
        }

    protected:
        double raw_pdf_(double x) const override final
        {
            const double size = arg(0);
            const double prob = arg(1);
            return exp(log_binom_coeff(size, x) + x * log(prob) + (size - x) * log1p(-prob));
        }
        double raw_cdf_(double x) const override final
        {
            const double size = arg(0);
            if (x > size)
                return 1.;
            return boost::math::ibeta(size - x, x + 1, 1. - arg(1), no_promote_policy{});
        }
        double raw_log_pdf_(double x) const override final
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
        double raw_rand_(rnd_gen &prng) const override final
        {
            return std::binomial_distribution<unsigned>(arg(0), arg(1))(prng);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Binomial_, DiscreteDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////  *** NegativeBinomial Distribution : class ***  ///////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class NegativeBinomial_ final : public Distribution
    {
    public:
        NegativeBinomial_()
            : Distribution({{"size", na::DOUBLE}, {"prob", na::DOUBLE}}, 0.) { discr_mode_ = Mode::Discrete::Defined; }
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(NegativeBinomial_)
        HOUSEPI_OVERRIDE_NAME_METHOD(negative_binomial)

        INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            const double size = arg(0);
            const double prob = arg(1);
            return exp(log_binom_coeff(x + size - 1., x) + x * log1p(-prob) + size * log(prob));
        }
        double raw_cdf_(double x) const override final
        {
            return boost::math::ibeta(arg(0), x + 1., arg(1));
        }
        double raw_log_pdf_(double x) const override final
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
        double raw_rand_(rnd_gen &prng) const override final
        {
            return std::negative_binomial_distribution<unsigned>(arg(0), arg(1))(prng);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(NegativeBinomial_, DiscreteDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Geometric Distribution : class ***  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Geometric_ final : public Distribution
    {
    public:
        Geometric_()
            : Distribution({{"prob", na::DOUBLE}}, 0.) { discr_mode_ = Mode::Discrete::Defined; }
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Geometric_)
        HOUSEPI_OVERRIDE_NAME_METHOD(geometric)

        INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            const double prob = arg(0);
            return prob * exp(x * log1p(-prob));
        }
        double raw_cdf_(double x) const override final
        {
            return 1. - exp((x + 1) * log1p(-arg(0)));
        }
        double raw_log_pdf_(double x) const override final
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
            p *= 1. - 64. * std::numeric_limits<double>::epsilon();
            return std::ceil(log1p(-p) / log1p(-arg(0)) - 1.);
        }
        double raw_rand_(rnd_gen &prng) const override final
        {
            return std::geometric_distribution<unsigned>(arg(0))(prng);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Geometric_, DiscreteDistrib)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Categorical Distribution : class ***  /////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Categorical_ final : public Distribution
    {
    public:
        Categorical_() : Distribution(), unif_gen(0., 1.) { discr_mode_ = Mode::Discrete::Defined; }
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Categorical_)
        HOUSEPI_OVERRIDE_NAME_METHOD(categorical)

        void init_discretization() override final {}
        void init_truncation() override final
        {
            const size_t cases = n_args();
            if (cases == 0)
                msg::stop(
                    "Invoked " + err_name() + " with no weights.\nProvide at least one weight.");

            if (weight_sum_() <= 0.)
                msg::stop(
                    "Invoked " + err_name() + " with weights summing to non-strictly positive value.");

            lower_.init_numeric(0.);
            upper_.init_numeric(static_cast<double>(cases - 1));
            below_mode_ = Mode::Interval::Closed;
            above_mode_ = Mode::Interval::Closed;
        }

    private:
        /// Live normalizing constant.
        double weight_sum_() const
        {
            double z = 0.;
            for (size_t i = 0; i < n_args(); ++i)
                z += arg(i);
            return z;
        }

    protected:
        Attribute &find_arg_(std::string_view key) override final { return find_or_append_arg_(key); }

        mutable std::uniform_real_distribution<double> unif_gen;

        double raw_pdf_(double x) const override final
        {
            return arg(static_cast<size_t>(x)) / weight_sum_();
        }
        double raw_log_pdf_(double x) const override final
        {
            return std::log(arg(static_cast<size_t>(x))) - std::log(weight_sum_());
        }
        double raw_cdf_(double x) const override final
        {
            const size_t k = static_cast<size_t>(x);
            double acc = 0.;
            for (size_t i = 0; i <= k; ++i)
                acc += arg(i);
            return acc / weight_sum_();
        }
        double raw_log_ratio_(double xo, double xn) const override final
        {
            return std::log(arg(static_cast<size_t>(xn))) - std::log(arg(static_cast<size_t>(xo)));
        }
        double raw_quantile_(double p) const override final
        {
            const size_t cases = n_args();
            const double target = p * weight_sum_();
            size_t i = 0;
            double acc = arg(0);
            while (i < cases && acc < target)
                acc += arg(++i);
            return static_cast<double>(i);
        }
        double raw_rand_(rnd_gen &prng) const override final
        {
            const size_t cases = n_args();
            const double target = unif_gen(prng) * weight_sum_();
            size_t i = 0;
            double acc = arg(0);
            while (i < cases && acc < target)
                acc += arg(++i);
            return static_cast<double>(i);
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Categorical_, DiscreteDistrib)

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** Ferretti20_ : class ***  ///////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    /// The infectivity profile is considered equal to the distribution
    /// of the Time from Onset of Symptoms to Transmission 't_ost':
    /// >	t_ost = <time infection infectee> - <time onset infector>
    /// Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
    /// where the pre-symptomatic values are rescaled as a function of the incubation time 'incub'.
    /// The arguments to provide are :
    ///   0) "shape";
    ///   1) "scale";
    ///   2) "c_loc" (custom location);
    ///   3) "mean_incub" (mean incubation period);
    ///   4) "incub_period" (individual incubation period);
    /// Note: symtoms can occur after infection (pre-symptomatic transmission),
    /// i. e. 't_ost' can be negative.
    /// Based on @cite ferretti_timing_2020
    class Ferretti20_ final : public Distribution
    {
    public:
        Ferretti20_()
            : Distribution({{"shape", 5.85}, {"scale", 1.85}, {"c_loc", -4.0}, {"mean_incub", 5.42}, {"incub_period", na::DOUBLE}}, 0.) {}
        HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Ferretti20_)
        HOUSEPI_OVERRIDE_NAME_METHOD(ferretti20)

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            const double a = arg(0); // 'shape'
            const double s = arg(1); // 'scale'
            const double m = arg(2); // 'c_loc'
            const double T = arg(3); // 'mean_incub'
            const double i = arg(4); // 'incub_period'
            double tOST = x - i;     // time from onset to transmission

            const double iT = i / T;
            const double norm_pdf = (s / a) * (1. + (iT - 1.) * pow(1 + exp(m / s), -a) - iT * pow(1. + exp((m + T) / s), -a));

            if (tOST < 0)
                tOST /= iT;

            const double emOSTs = exp((m - tOST) / s);
            return emOSTs * pow(1. + emOSTs, -a - 1.) / norm_pdf;
        }

        double raw_cdf_(double x) const override final
        {
            const double a = arg(0);   // 'shape'
            const double s = arg(1);   // 'scale'
            const double m = arg(2);   // 'c_loc'
            const double T = arg(3);   // 'mean_incub'
            const double i = arg(4);   // 'incub_period'
            const double tOST = x - i; // time from onset to transmission

            const double iT = i / T;
            const double pema = pow(1. + exp(m / s), -a);
            const double pemTSa = pow(1. + exp((m + T) / s), -a);
            const double norm_cdf = 1. + (iT - 1.) * pema - iT * pemTSa;

            if (tOST < 0) // Integrate from -incub_period to tOST
                return iT * (pow(1. + exp((m - tOST / iT) / s), -a) - pemTSa) / norm_cdf;

            // Integrate from -incub to 0
            const double F0 = iT * (pema - pemTSa);
            // Integrate from 0 to t_ost
            return (F0 + pow(1. + exp((m - tOST) / s), -a) - pema) / norm_cdf;
        }
    };
    HOUSEPI_REGISTER_FUNCTION(Ferretti20_, InfecProfile)
}
