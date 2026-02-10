/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
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

#include "tools/string.h"

#include "distribution.h"
#include "registry.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Distribution - reshaping  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Select functional method in case of truncation.
void Distribution::helper_init_distr_truncation_()
{
    trunc_mode_ = Trunc::Active;

    if (lower_.is_empty() && upper_.is_empty()) // No Truncation
    {
        trunc_ = &Distribution::do_trunc_dummy_;
        shift_ = &Distribution::do_shift_dummy_;
        trunc_mode_ = Trunc::None;
    }
    else if (!lower_.is_empty() && upper_.is_empty()) // Lower Range Truncation
    {
        switch (discr_mode_)
        {
        case Mode::Discrete::Defined: // discrete defined
            trunc_ = &Distribution::do_trunc_min_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;

            break;
        case Mode::Discrete::Center: // centered discretization
            trunc_ = &Distribution::do_trunc_min_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            break;
        case Mode::Discrete::Left: // left discretization
            trunc_ = &Distribution::do_trunc_min_discr_left_;
            shift_ = &Distribution::do_shift_min_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_min_;
            shift_ = &Distribution::do_shift_min_;
        }
    }
    else if (lower_.is_empty() && !upper_.is_empty()) // Upper Range Truncation
    {
        shift_ = &Distribution::do_shift_dummy_;
        switch (discr_mode_)
        {
        case Mode::Discrete::Defined: // discrete defined
            trunc_ = &Distribution::do_trunc_max_discr_defined_;
            break;
        case Mode::Discrete::Center: // centered discretization
            trunc_ = &Distribution::do_trunc_max_discr_center_;
            break;
        case Mode::Discrete::Left: // left discretization
            trunc_ = &Distribution::do_trunc_max_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_max_;
        }
    }
    else // Both Sides Truncation
    {
        switch (discr_mode_)
        {
        case Mode::Discrete::Defined: // discrete defined
            trunc_ = &Distribution::do_trunc_both_sides_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            break;
        case Mode::Discrete::Center: // centered discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            break;
        case Mode::Discrete::Left: // left discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_left_;
            shift_ = &Distribution::do_shift_min_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_both_sides_;
            shift_ = &Distribution::do_shift_min_;
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Distribution Registry : getter  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<Function> DistributionRegistry::create(std::string_view opt, bool is_prim) const
{
    auto it = factory_map_.find(opt);
    if (it == factory_map_.end())
        return nullptr;

    auto output = it->second();
    output->set_eval_as_distr_call();
    if (is_prim)
        output->set_eval_as_distr_primitive();

    return output;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Distribution Registry : setter  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<Distribution> compose_distr_registries(std::string_view opt,
                                                       const std::string &err_type,
                                                       const std::vector<std::string> &registries,
                                                       bool is_cumul)
{
    if (registries.empty())
        throw std::invalid_argument(
            "No registries were provided for '" + err_type + "'.");

    if (hsp::na::is(opt))
        return nullptr;

    std::vector<const DistributionRegistry *> registry_ptrs = {};
    for (const std::string &r : registries)
        if (r == "continuos")
            registry_ptrs.push_back(&ContinuosDistrRegistry::instance());
        else if (r == "discrete")
            registry_ptrs.push_back(&DiscreteDistrRegistry::instance());
        else if (r == "infecprofile")
            registry_ptrs.push_back(&InfecProfileRegistry::instance());
        else
            throw std::invalid_argument(
                "Unknown " + err_type + " registry '" + r + "'.");

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

    throw std::invalid_argument(
        "Unknown " + err_type + " '" + std::string(opt) + "'." +
        hsp::str::find_closest(opt, options));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  *** Distribution : implementations ***  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
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
#define INIT_DISCR_LEFT_MIN_ZERO_CLOSE()       \
    void init_discretization() override final  \
    {                                          \
        discr_mode_ = Mode::Discrete::Left;    \
    }                                          \
    void init_truncation() override final      \
    {                                          \
        if (helper_init_range_min_(0., false)) \
            helper_init_distr_truncation_();   \
    }

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
        std::string_view name() const noexcept override final { return "uniform"; }

        /// UPDATE: should understand when discrete defined...
        void init_discretization() override final { discr_mode_ = Mode::Discrete::Center; }

    protected:
        mutable std::uniform_real_distribution<double> unif_gen;

        double raw_pdf_(double) const override final { return 1. / (upper_.value() - lower_.value()); }
        double raw_cdf_(double x) const override final { return (x - lower_.value()) / (upper_.value() - lower_.value()); }
        double raw_logpdf_(double) const override final { return 0.; }
        double raw_logratio_(double, double) const override final { return 0.; }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            const double r = unif_gen(prng);
            return lower_.value() * (1. - r) + upper_.value() * r;
        }
        double raw_quantile_(double p) const override final { return lower_.value() * (1. - p) + upper_.value() * p; }
    };
    REGISTER_HOUSEPI_FUNCTION(Uniform_, ContinuosDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Normal Distribution : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Normal_ final : public Distribution
    {
    public:
        Normal_()
            : Distribution({{"mean", 0.}, {"sd", 1.}}),
              norm_gen(0., 1.) {}
        std::string_view name() const noexcept override final { return "normal"; }

        void init_discretization() override final
        {
            discr_mode_ = lower_.is_empty() ? Mode::Discrete::Center : Mode::Discrete::Left;
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
        double raw_logpdf_(double x) const override final
        {
            const double sd = arg(1);
            const double z = (x - arg(0)) / sd;
            return -0.5 * z * z - M_LN_SQRT_2PI - log(sd);
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            const double sd = arg(1);
            return (0.5 * xo * xo + (xn - xo) * arg(0) - 0.5 * xn * xn) / (sd * sd);
        }
        double raw_quantile_(double p) const override final
        {
            return arg(0) + arg(1) * boost::math::erfc_inv(2. - 2. * p, no_promote_policy{}) * M_SQRT2;
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return arg(0) + arg(1) * norm_gen(prng);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Normal_, ContinuosDistrRegistry)

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
        std::string_view name() const noexcept override final { return "lognormal"; }

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

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
        double raw_logpdf_(double x) const override final
        {
            const double sdlog = arg(1);
            const double z = (log(x) - arg(0)) / sdlog;
            return -0.5 * z * z - log(x * sdlog) - M_LN_SQRT_2PI;
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            const double lxo = log(xo);
            const double lxn = log(xn);
            const double meanlog = arg(0);
            const double sdlog = arg(1);
            return 0.5 * ((lxo - meanlog) * lxo - (lxn - meanlog) * lxn) / (sdlog * sdlog) + lxo - lxn;
        }
        double raw_quantile_(double p) const override final
        {
            return exp(arg(0) + arg(1) * boost::math::erfc_inv(2. - 2. * p, no_promote_policy{}) * M_SQRT2);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return exp(arg(0) + arg(1) * norm_gen(prng));
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Lognormal_, ContinuosDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** Beta Distribution : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Beta_ final : public Distribution
    {
    public:
        Beta_()
            : Distribution({{"shape1", hsp::na::DOUBLE}, {"shape2", hsp::na::DOUBLE}},
                           0., 1.)
        {
            below_mode_ = Mode::Interval::Closed;
            above_mode_ = Mode::Interval::Closed;
        }
        std::string_view name() const noexcept override final { return "beta"; }

        void init_discretization() override final
        {
            throw std::invalid_argument(
                "Invoked invalid option discrete for " + err_name() + ".");
        }
        void init_truncation() override final
        {
            const double sx = lower_.value();
            const double dx = upper_.value();
            if (sx < 0)
                throw std::invalid_argument(
                    "Invoked invalid min value " + hsp::str::format(sx) + " for " + err_name() + ".");
            if (dx > 1)
                throw std::invalid_argument(
                    "Invoked invalid max value " + hsp::str::format(dx) + " for " + err_name() + ".");
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
        double raw_logpdf_(double x) const override final
        {
            const double alpha = arg(0);
            const double beta = arg(1);
            return (alpha - 1.) * log(x) + (beta - 1.) * log1p(-x) - lbeta(alpha, beta);
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            return (arg(0) - 1.) * log(xn / xo) + (arg(1) - 1.) * log((1. - xn) / (1. - xo));
        }
        double raw_quantile_(double p) const override final
        {
            return boost::math::ibeta_inv(arg(0), arg(1), p, no_promote_policy{});
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return boost::random::beta_distribution<double>(arg(0), arg(1))(prng);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Beta_, ContinuosDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** Gamma Distribution : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Gamma_ final : public Distribution
    {
    public:
        Gamma_()
            : Distribution({{"shape", hsp::na::DOUBLE}, {"scale", 1.}}, 0.) {}
        std::string_view name() const noexcept override final { return "gamma"; }

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            if (x != 0.)
                return exp(raw_logpdf_(x));

            const double shape = arg(0);
            if (shape < 1.)
                return std::numeric_limits<double>::infinity();
            else if (shape == 1.)
                return exp(-log(arg(1)) - lgamma(1));
            else
                return 0.;
        }
        double raw_cdf_(double x) const override final
        {
            return boost::math::gamma_p(
                arg(0), x / arg(1), no_promote_policy{});
        }
        double raw_logpdf_(double x) const override final
        {
            const double shape = arg(0);
            const double t = x / arg(1);
            return shape * log(t) - log(x) - lgamma(shape) - t;
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            return (arg(0) - 1.) * log(xn / xo) + (xo - xn) / arg(1);
        }
        double raw_quantile_(double p) const override final
        {
            return arg(1) * boost::math::gamma_p_inv(
                                arg(0), p, no_promote_policy{});
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::gamma_distribution<double>(arg(0), arg(1))(prng);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Gamma_, ContinuosDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////  *** Gamma (Mean-Varance) Distribution : class ***  /////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Gamma (Mean-Variance) Distribution.
    /// UPDATE: this should be collapsed with the previous
    class GammaMeanVar_ final : public Distribution
    {
    public:
        GammaMeanVar_()
            : Distribution({{"mean", hsp::na::DOUBLE}, {"var", hsp::na::DOUBLE}}, 0.) {}
        std::string_view name() const noexcept override final { return "gamma_meanvar"; }

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            if (x != 0.)
                return exp(raw_logpdf_(x));

            const double m = arg(0);
            const double v = arg(1);
            const double shape = m * m / v;
            if (shape < 1.)
                return std::numeric_limits<double>::infinity();
            else if (shape == 1.)
                return exp(-log(v / m) - lgamma(1));
            else
                return 0.;
        }
        double raw_cdf_(double x) const override final
        {
            const double m = arg(0);
            const double v = arg(1);
            return boost::math::gamma_p(
                m * m / v, x * m / v, no_promote_policy{});
        }
        double raw_logpdf_(double x) const override final
        {
            const double m = arg(0);
            const double v = arg(1);
            const double shape = m * m / v;
            const double t = x * m / v;
            return shape * log(t) - log(x) - lgamma(shape) - t;
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            const double m = arg(0);
            const double v = arg(1);
            return (m * m / v - 1.) * log(xn / xo) + (xo - xn) * m / v;
        }
        double raw_quantile_(double p) const override final
        {
            const double m = arg(0);
            const double v = arg(1);
            return v / m * boost::math::gamma_p_inv(m * m / v, p, no_promote_policy{});
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            const double m = arg(0);
            const double v = arg(1);
            return std::gamma_distribution<double>(m * m / v, v / m)(prng);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(GammaMeanVar_, ContinuosDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Rayleigh Distribution : class ***  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Rayleigh_ final : public Distribution
    {
    public:
        Rayleigh_()
            : Distribution({{"scale", 1.}}, 0.),
              unif_gen(0., 1.) {}
        std::string_view name() const noexcept override final { return "rayleigh"; }

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
            const double scale = arg(0);
            return 1. - exp(-0.5 * x * x / (scale * scale));
        }
        double raw_logpdf_(double x) const override final
        {
            const double scale = arg(0);
            const double scale2 = scale * scale;
            return log(x / scale2) - 0.5 * x * x / scale2;
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            const double scale = arg(0);
            return log(xn / xo) + 0.5 * (xo * xo - xn * xn) / (scale * scale);
        }
        double raw_quantile_(double p) const override final
        {
            return arg(0) * sqrt(-2. * log1p(-p));
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return arg(0) * sqrt(-2 * log(unif_gen(prng)));
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Rayleigh_, ContinuosDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Weibull Distribution : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Weibull_ final : public Distribution
    {
    public:
        Weibull_()
            : Distribution({{"shape", hsp::na::DOUBLE}, {"scale", 1.}}, 0.),
              unif_gen(0., 1.) {}
        std::string_view name() const noexcept override final { return "weibull"; }

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
            return 1. - exp(-pow(x / arg(1), arg(0)));
        }
        double raw_logpdf_(double x) const override final
        {
            const double shape = arg(0);
            const double scale = arg(1);
            const double z = x / scale;
            return log(shape / scale) + (shape - 1.) * log(z) - pow(z, shape);
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            const double shape = arg(0);
            const double scale = arg(1);
            return (shape - 1.) * log(xn / xo) + pow(xo / scale, shape) - pow(xn / scale, shape);
        }
        double raw_quantile_(double p) const override final
        {
            return arg(1) * pow(-log1p(-p), 1. / arg(0));
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return arg(1) * pow(-log(unif_gen(prng)), 1.0 / arg(0));
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Weibull_, ContinuosDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  *** Exponential Distribution : class ***  //////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Exponential_ final : public Distribution
    {
    public:
        Exponential_()
            : Distribution({{"rate", 1.}}, 0.),
              unif_gen(0., 1.) {}
        std::string_view name() const noexcept override final { return "exponential"; }

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
            return 1. - exp(-arg(0) * x);
        }
        double raw_logpdf_(double x) const override final
        {
            const double rate = arg(0);
            return log(rate) - rate * x;
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            return arg(0) * (xo - xn);
        }
        double raw_quantile_(double p) const override final
        {
            return -log1p(-p) / arg(0);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return -log(unif_gen(prng)) / arg(0);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Exponential_, ContinuosDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** Poisson Distribution : class ***  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Poisson_ final : public Distribution
    {
    public:
        Poisson_()
            : Distribution({{"lambda", hsp::na::DOUBLE}}, 0.) { discr_mode_ = Mode::Discrete::Defined; }
        std::string_view name() const noexcept override final { return "poisson"; }

        INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            return exp(raw_logpdf_(x));
        }
        double raw_cdf_(double x) const override final
        {
            return 1. - boost::math::gamma_p(x + 1., arg(0), no_promote_policy{});
        }
        double raw_logpdf_(double x) const override final
        {
            const double lambda = arg(0);
            return x * log(lambda) - lambda - lgamma(x + 1.);
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            return (xn - xo) * log(arg(0)) + lgamma(xo + 1.) - lgamma(xn + 1.);
        }
        double raw_quantile_(double p) const override final
        {
            return boost::math::quantile(boost::math::poisson_distribution<double>(arg(0)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::poisson_distribution<unsigned>(arg(0))(prng);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Poisson_, DiscreteDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** Binomial Distribution : class ***  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Binomial_ final : public Distribution
    {
    public:
        Binomial_()
            : Distribution({{"size", hsp::na::DOUBLE}, {"prob", hsp::na::DOUBLE}}, 0.)
        {
            discr_mode_ = Mode::Discrete::Defined;
            upper_.init_pointer(arguments_[0].pointer());
            below_mode_ = Mode::Interval::Closed;
            above_mode_ = Mode::Interval::Closed;
        }
        std::string_view name() const noexcept override final { return "binomial"; }

        void init_discretization() override final {}
        void init_truncation() override final
        {
            const double sx = lower_.value();
            if (sx < 0.)
                throw std::invalid_argument(
                    "Invoked invalid min value " + hsp::str::format(sx) +
                    " for " + err_name() + ".");
            if (sx > 0. || arguments_[0].pointer() != upper().pointer())
            {
                const double dx = lower_.value();
                if (dx > arg(0))
                    throw std::invalid_argument(
                        "Invoked invalid max value " + hsp::str::format(dx) +
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
        double raw_logpdf_(double x) const override final
        {
            const double size = arg(0);
            const double prob = arg(1);
            return log_binom_coeff(size, x) + x * log(prob) + (size - x) * log1p(-prob);
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            const double prob = arg(1);
            return logratio_binom_coeff(arg(0), xo, xn) + (xn - xo) * log(prob / (1. - prob));
        }
        double raw_quantile_(double p) const override final
        {
            return boost::math::quantile(boost::math::binomial_distribution<double>(arg(0), arg(1)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::binomial_distribution<unsigned>(arg(0), arg(1))(prng);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Binomial_, DiscreteDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////  *** NegativeBinomial Distribution : class ***  ///////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class NegativeBinomial_ final : public Distribution
    {
    public:
        NegativeBinomial_()
            : Distribution({{"size", hsp::na::DOUBLE}, {"prob", hsp::na::DOUBLE}}, 0.) { discr_mode_ = Mode::Discrete::Defined; }
        std::string_view name() const noexcept override final { return "negative_binomial"; }

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
        double raw_logpdf_(double x) const override final
        {
            const double size = arg(0);
            const double prob = arg(1);
            return log_binom_coeff(x + size - 1., x) + x * log1p(-prob) + size * log(prob);
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            const double size = arg(0);
            const double prob = arg(1);
            return logratio_binom_coeff(size, xo, xn) + (xn - xo) * log1p(-prob);
        }
        double raw_quantile_(double p) const override final
        {
            return boost::math::quantile(boost::math::negative_binomial_distribution<double>(arg(0), arg(1)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::negative_binomial_distribution<unsigned>(arg(0), arg(1))(prng);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(NegativeBinomial_, DiscreteDistrRegistry)

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Geometric Distribution : class ***  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Geometric_ final : public Distribution
    {
    public:
        Geometric_()
            : Distribution({{"prob", hsp::na::DOUBLE}}, 0.) { discr_mode_ = Mode::Discrete::Defined; }
        std::string_view name() const noexcept override final { return "geometric"; }

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
        double raw_logpdf_(double x) const override final
        {
            const double prob = arg(0);
            return log(prob) + x * log1p(-prob);
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            return (xn - xo) * log1p(-arg(0));
        }
        double raw_quantile_(double p) const override final
        {
            return std::floor(log1p(-p) / log1p(-arg(0)));
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::geometric_distribution<unsigned>(arg(0))(prng);
        }
    };
    REGISTER_HOUSEPI_FUNCTION(Geometric_, DiscreteDistrRegistry)

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** InfecProfile_Ferretti20 : class ***  ///////////////////////////
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
    ///
    /// Ref:
    /// ----
    /// Ferretti L, Ledda A, Wymant C, Zhao L, Ledda V, Abeler-Dorner L, et al.
    /// The timing of COVID-19 transmission.
    /// medRxiv; 2020. p. 2020.09.04.20188516. doi:10.1101/2020.09.04.20188516
    class InfecProfile_Ferretti20 final : public Distribution
    {
    public:
        InfecProfile_Ferretti20()
            : Distribution({{"shape", 5.85}, {"scale", 1.85}, {"c_loc", -4.0}, {"mean_incub", 5.42}, {"incub_period", 5.42}}, 0.) {}
        std::string_view name() const noexcept override final { return "ferretti20"; }

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

            const double norm_pdf = (s / a) * (1. + (i / T - 1.) * pow(1 + exp(m / s), -a) - (i / T) * pow(1. + exp((m + T) / s), -a));

            if (tOST < 0)
                tOST *= (T / i);

            return exp((m - tOST) / s) * pow(1. + exp((m - tOST) / s), -a - 1.) / norm_pdf;
        }

        double raw_cdf_(double x) const override final
        {
            const double a = arg(0);   // 'shape'
            const double s = arg(1);   // 'scale'
            const double m = arg(2);   // 'c_loc'
            const double T = arg(3);   // 'mean_incub'
            const double i = arg(4);   // 'incub_period'
            const double tOST = x - i; // time from onset to transmission

            const double norm_cdf = (1. + (i / T - 1.) * pow(1 + exp(m / s), -a) - (i / T) * (pow(1. + exp((m + T) / s), -a)));

            if (tOST < 0) // Integrate from -incub_period to tOST
                return (i / T) * (pow(1. + exp((m - tOST * T / i) / s), -a) - pow(1. + exp((m + T) / s), -a)) / norm_cdf;

            // Integrate from -incub to 0
            const double F0 = (i / T) * (pow(1. + exp(m / s), -a) - pow(1. + exp((m + T) / s), -a));
            // Integrate from 0 to t_ost
            return (F0 + pow(1. + exp((m - tOST) / s), -a) - pow(1. + exp(m / s), -a)) / norm_cdf;
        }
    };
    REGISTER_HOUSEPI_FUNCTION(InfecProfile_Ferretti20, InfecProfileRegistry)
}
