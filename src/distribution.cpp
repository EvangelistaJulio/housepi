/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/random/beta_distribution.hpp>

#include "tools/string.h"

#include "distribution.h"
#include "registry.h"

namespace
{
    using std::string;
    using std::vector;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Distribution - reshaping  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Select functional for discretization of random sampling.
void Distribution::init_discr_helper_(unsigned short discr_type)
{
    Function::init_discr_helper_(discr_type);

    if (discr_type == hsp::discr::CENTER) // centered discretization
    {
        rand_ = &Distribution::do_discr_center_rand_;
        if (!loc_arg_.is_empty())
            rand_ = &Distribution::do_discr_center_rand_loc_;
    }
    else if (discr_type == hsp::discr::LEFT) // left discretization
    {
        rand_ = &Distribution::do_discr_left_rand_;
        if (!loc_arg_.is_empty())
            rand_ = &Distribution::do_discr_left_rand_loc_;
    }
    else // discrete defined or no discretization
    {
        rand_ = &Distribution::raw_rand_;
        if (!loc_arg_.is_empty())
            rand_ = &Distribution::do_rand_loc_;
    }
}

/// Select functional method in case of truncation.
void Distribution::init_trunc_helper_()
{
    const unsigned short discr_type = which_discrete_case_();

    if (lower_arg_.is_empty() && upper_arg_.is_empty()) // No Truncation
    {
        trunc_ = &Distribution::do_trunc_dummy_;
        logtrunc_ = &Distribution::do_logtrunc_dummy_;
        shift_ = &Distribution::do_shift_dummy_;
    }
    else if (!lower_arg_.is_empty() && upper_arg_.is_empty()) // Lower Range Truncation
    {
        switch (discr_type)
        {
        case hsp::discr::DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_min_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case hsp::discr::CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_min_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case hsp::discr::LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_min_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_left_;
            shift_ = &Distribution::do_shift_min_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_min_;
            logtrunc_ = &Distribution::do_logtrunc_min_;
            shift_ = &Distribution::do_shift_min_;
            rand_ = &Distribution::do_trunc_rand_;
        }
    }
    else if (lower_arg_.is_empty() && !upper_arg_.is_empty()) // Upper Range Truncation
    {
        shift_ = &Distribution::do_shift_dummy_;
        switch (discr_type)
        {
        case hsp::discr::DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_max_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case hsp::discr::CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_max_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case hsp::discr::LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_max_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_max_;
            logtrunc_ = &Distribution::do_logtrunc_max_;
            rand_ = &Distribution::do_trunc_rand_;
        }
    }
    else // Both Sides Truncation
    {
        switch (discr_type)
        {
        case hsp::discr::DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_both_sides_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case hsp::discr::CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case hsp::discr::LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_left_;
            shift_ = &Distribution::do_shift_min_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_both_sides_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_;
            shift_ = &Distribution::do_shift_min_;
            rand_ = &Distribution::do_trunc_rand_;
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
                                                       const string &err_type,
                                                       const vector<std::string> &registries,
                                                       bool is_cumul)
{
    if (registries.empty())
        throw std::invalid_argument(
            "No registries were provided for '" + err_type + "'.");

    if (hsp::na::is(opt))
        return nullptr;

    vector<const DistributionRegistry *> registry_ptrs = {};
    for (const auto &r : registries)
        if (r == "continuos")
            registry_ptrs.push_back(&ContinuosDistrRegistry::instance());
        else if (r == "discrete")
            registry_ptrs.push_back(&DiscreteDistrRegistry::instance());
        else if (r == "infecprofile")
            registry_ptrs.push_back(&InfecProfileRegistry::instance());
        else
            throw std::invalid_argument("Unknown " + err_type + " registry '" + r + "'.");

    vector<string> options = {};
    for (auto &r : registry_ptrs)
    {
        if (r->count(opt))
            return std::unique_ptr<Distribution>(
                static_cast<Distribution *>(r->create(opt, is_cumul).release()));
        const auto &y = r->get_options();
        options.reserve(options.size() + y.size());
        options.insert(options.end(), y.begin(), y.end());
    }

    throw std::invalid_argument(
        "Unknown " + err_type + " '" + string(opt) + "'." + hsp::str::find_closest(opt, options));
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
    constexpr double M_1_SQRT_2PI = 0.398942280401432677939946059934; /// 1/sqrt(2pi
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
#define INIT_DISCR_LEFT_MIN_ZERO_CLOSE()      \
    void init_discretization() override final \
    {                                         \
        init_discr_helper_(hsp::discr::LEFT); \
    }                                         \
    void init_truncation() override final     \
    {                                         \
        if (init_range_min_(0., false))       \
            init_trunc_helper_();             \
    }

// Macro for discrete defined when min = 0 (close).
#define INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()      \
    void init_discretization() override final {} \
    void init_truncation() override final        \
    {                                            \
        if (init_range_min_(0., false))          \
            init_trunc_helper_();                \
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  *** Uniform Distribution : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Uniform_ : public Distribution
    {
    public:
        Uniform_()
            : Distribution("uniform",
                           {{hsp::params::MIN, 0.}, {hsp::params::MAX, 1.}},
                           hsp::na::DOUBLE, false, hsp::na::DOUBLE, false),
              unif_gen(arg(0), arg(1)) { init_uniform_range_(); }

        void init_discretization() override final
        {
            // UPDATE: should understand when discrete defined...
            init_discr_helper_(hsp::discr::CENTER);
        }
        void init_truncation() override final
        {
            init_uniform_range_();
            init_range_helper_(false, false);
            // NOTE: truncation for uniform distribution is inside raw
        }

    protected:
        mutable std::uniform_real_distribution<double> unif_gen;

        double raw_pdf_(double) const override final { return 1. / (arg(1) - arg(0)); }
        double raw_cdf_(double x) const override final { return (x - arg(0)) / (arg(1) - arg(0)); }
        double raw_logpdf_(double) const override final { return 0.; }
        double raw_logratio_(double, double) const override final { return 0.; }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            const double r = unif_gen(prng);
            return arg(0) * (1. - r) + arg(1) * r;
        }
        double raw_quantile_(double p) const override final { return arg(0) * (1. - p) + arg(1) * p; }

    private:
        /// Ensure uniform range (min, max) match argument (min, max).
        void init_uniform_range_()
        {
            // impose range minimum equal to argument minimum
            const auto &min_arg = arguments_[0];
            if (min_arg.is_function())
                lower_arg_.transfer(min_arg);
            else
                lower_arg_.transfer(min_arg);

            // impose range maximum equal to argument maximum
            const auto &max_arg = arguments_[1];
            if (max_arg.is_function())
                upper_arg_.transfer(max_arg);
            else
                upper_arg_.transfer(max_arg);
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Normal Distribution : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Normal_ : public Distribution
    {
    public:
        Normal_()
            : Distribution("normal",
                           {{"mean", 0.}, {"sd", 1.}},
                           hsp::na::DOUBLE, false, hsp::na::DOUBLE, false),
              norm_gen(0., 1.) {}

        void init_discretization() override final
        {
            if (lower_arg_.is_empty())
                init_discr_helper_(hsp::discr::CENTER);
            else
                init_discr_helper_(hsp::discr::LEFT);
        }
        void init_truncation() override final
        {
            init_range_helper_(false, false);
            init_trunc_helper_();
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
            static const boost::math::normal stdnorm(0., 1.);
            return arg(0) + arg(1) * boost::math::quantile(stdnorm, p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return arg(0) + arg(1) * norm_gen(prng);
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  *** Lognormal Distribution : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Lognormal Distribution derived class.
    class Lognormal_ : public Distribution
    {
    public:
        Lognormal_()
            : Distribution("lognormal",
                           {{"meanlog", 0.}, {"sdlog", 1.}},
                           0., true, hsp::na::DOUBLE, false),
              norm_gen(0., 1.) {}

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
            static const boost::math::normal stdnorm(0., 1.);
            return exp(arg(0) + arg(1) * boost::math::quantile(stdnorm, p));
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return exp(arg(0) + arg(1) * norm_gen(prng));
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** Beta Distribution : class ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Beta_ : public Distribution
    {
    public:
        Beta_()
            : Distribution("beta",
                           {{"shape1", hsp::na::DOUBLE}, {"shape2", hsp::na::DOUBLE}},
                           0., true, 1., true) {}

        void init_discretization() override final
        {
            throw std::invalid_argument(
                "Invoked invalid option discrete for " + err_name() + ".");
        }
        void init_truncation() override final
        {
            if (lower() < 0)
                throw std::invalid_argument(
                    "Invoked invalid min value " + hsp::str::format(lower()) + " for " + err_name() + ".");
            if (upper() > 1)
                throw std::invalid_argument(
                    "Invoked invalid max value " + hsp::str::format(upper()) + " for " + err_name() + ".");
            if (lower() > 0 || upper() < 1)
            {
                init_range_helper_(lower() == 0, upper() == 1);
                init_trunc_helper_();
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
            return boost::math::ibeta(arg(0), arg(1), x);
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
            return boost::math::quantile(boost::math::beta_distribution<>(arg(0), arg(1)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return boost::random::beta_distribution<double>(arg(0), arg(1))(prng);
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////  *** Gamma Distribution : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Gamma_ : public Distribution
    {
    public:
        Gamma_()
            : Distribution("gamma",
                           {{"shape", hsp::na::DOUBLE}, {"scale", 1.}},
                           0., false, hsp::na::DOUBLE, false) {}

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            const double shape = arg(0);
            const double scale = arg(1);
            return pow(x, shape - 1.) * exp(-shape * log(scale) - lgamma(shape) - x / scale);
        }
        double raw_cdf_(double x) const override final
        {
            return boost::math::gamma_p(arg(0), x / arg(1));
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
            return boost::math::quantile(boost::math::gamma_distribution<>(arg(0), arg(1)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::gamma_distribution<double>(arg(0), arg(1))(prng);
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////  *** Gamma (Mean-Varance) Distribution : class ***  /////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Gamma (Mean-Variance) Distribution.
    /// UPDATE: this should be collapsed with the previous
    class GammaMeanVar_ : public Distribution
    {
    public:
        GammaMeanVar_()
            : Distribution("gamma_meanvar",
                           {{"mean", hsp::na::DOUBLE}, {"var", hsp::na::DOUBLE}},
                           0., false, hsp::na::DOUBLE, false) {}

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            const double m = arg(0);
            const double v = arg(1);
            const double shape = m * m / v;
            const double scale = v / m;
            return pow(x, shape - 1.) * exp(-shape * log(scale) - lgamma(shape) - x / scale);
        }
        double raw_cdf_(double x) const override final
        {
            const double m = arg(0);
            const double v = arg(1);
            return boost::math::gamma_p(m * m / v, x * m / v);
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
            return boost::math::quantile(boost::math::gamma_distribution<>(m * m / v, v / m), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            const double m = arg(0);
            const double v = arg(1);
            return std::gamma_distribution<double>(m * m / v, v / m)(prng);
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Rayleigh Distribution : class ***  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Rayleigh_ : public Distribution
    {
    public:
        Rayleigh_()
            : Distribution("rayleigh",
                           {{"scale", 1.}},
                           0., false, hsp::na::DOUBLE, false),
              unif_gen(0., 1.) {}

        INIT_DISCR_LEFT_MIN_ZERO_CLOSE()

    protected:
        mutable std::uniform_real_distribution<double> unif_gen;

        double raw_pdf_(double x) const override final
        {
            const double scale = arg(0);
            const double scale2 = scale * scale;
            return x * exp(-.5 * x * x / scale2) / scale2;
        }
        double raw_cdf_(double x) const override final
        {
            const double scale = arg(0);
            return 1. - exp(-.5 * x * x / (scale * scale));
        }
        double raw_logpdf_(double x) const override final
        {
            const double scale = arg(0);
            const double scale2 = scale * scale;
            return log(x / scale2) - .5 * x * x / scale2;
        }
        double raw_logratio_(double xo, double xn) const override final
        {
            const double scale = arg(0);
            return log(xn / xo) + .5 * (xo * xo - xn * xn) / (scale * scale);
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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Weibull Distribution : class ***  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Weibull_ : public Distribution
    {
    public:
        Weibull_()
            : Distribution("weibull",
                           {{"shape", hsp::na::DOUBLE}, {"scale", 1.}},
                           0., false, hsp::na::DOUBLE, false),
              unif_gen(0., 1.) {}

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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////  *** Exponential Distribution : class ***  //////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Exponential_ : public Distribution
    {
    public:
        Exponential_()
            : Distribution("exponential",
                           {{"rate", 1.}},
                           0., false, hsp::na::DOUBLE, false),
              unif_gen(0., 1.) {}

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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** Poisson Distribution : class ***  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Poisson_ : public Distribution
    {
    public:
        Poisson_()
            : Distribution("poisson",
                           {{"lambda", hsp::na::DOUBLE}},
                           0., false, hsp::na::DOUBLE, false) { init_discr_helper_(hsp::discr::DEFINED); }

        INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            const double lambda = arg(0);
            return exp(x * log(lambda) - lambda - lgamma(x + 1.));
        }
        double raw_cdf_(double x) const override final
        {
            return 1. - boost::math::gamma_p(x + 1., arg(0));
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
            return boost::math::quantile(boost::math::poisson_distribution<>(arg(0)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::poisson_distribution<unsigned>(arg(0))(prng);
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  *** Binomial Distribution : class ***  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Binomial_ : public Distribution
    {
    public:
        Binomial_()
            : Distribution("binomial",
                           {{"size", hsp::na::DOUBLE}, {"prob", hsp::na::DOUBLE}},
                           0., false, hsp::na::DOUBLE, false)
        {
            init_discr_helper_(hsp::discr::DEFINED);
            upper_arg_.init_pointer(get_arg(0).get_ptr());
            init_range_helper_(false, false);
        }

        void init_discretization() override final {}
        void init_truncation() override final
        {
            if (lower() < 0.)
                throw std::invalid_argument(
                    "Invoked invalid min value " + hsp::str::format(lower_arg_.get()) +
                    " for " + err_name() + ".");
            if (lower() > 0. || get_arg(0).get_ptr() != get_upper().get_ptr())
            {
                if (upper() > arg(0))
                    throw std::invalid_argument(
                        "Invoked invalid max value " + hsp::str::format(upper_arg_.get()) +
                        " for " + err_name() + " with size " + std::to_string(arg(0)) + ".");
                init_range_helper_(false, false);
                init_trunc_helper_();
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
            return boost::math::ibeta(size - x, x + 1, 1. - arg(1));
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
            return boost::math::quantile(boost::math::binomial_distribution<>(arg(0), arg(1)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::binomial_distribution<unsigned>(arg(0), arg(1))(prng);
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////  *** NegativeBinomial Distribution : class ***  ///////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class NegativeBinomial_ : public Distribution
    {
    public:
        NegativeBinomial_()
            : Distribution("negative_binomial",
                           {{"size", hsp::na::DOUBLE}, {"prob", hsp::na::DOUBLE}},
                           0., false, hsp::na::DOUBLE, false) { init_discr_helper_(hsp::discr::DEFINED); }

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
            return boost::math::quantile(boost::math::negative_binomial_distribution<>(arg(0), arg(1)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::negative_binomial_distribution<unsigned>(arg(0), arg(1))(prng);
        }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  *** Geometric Distribution : class ***  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    class Geometric_ : public Distribution
    {
    public:
        Geometric_()
            : Distribution("geometric",
                           {{"prob", hsp::na::DOUBLE}},
                           0., false, hsp::na::DOUBLE, false) { init_discr_helper_(hsp::discr::DEFINED); }

        INIT_DISCR_DEFINED_MIN_ZERO_CLOSE()

    protected:
        double raw_pdf_(double x) const override final
        {
            const double prob = arg(0);
            return prob * pow(1. - prob, x);
        }
        double raw_cdf_(double x) const override final
        {
            return 1. - pow(1. - arg(0), x + 1);
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
            return boost::math::quantile(boost::math::geometric_distribution<>(arg(0)), p);
        }
        double raw_rand_(hsp::rnd_gen &prng) const override final
        {
            return std::geometric_distribution<unsigned>(arg(0))(prng);
        }
    };

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
    class InfecProfile_Ferretti20 : public Distribution
    {
    public:
        InfecProfile_Ferretti20()
            : Distribution("ferretti20",
                           {{"shape", 5.85}, {"scale", 1.85}, {"c_loc", -4.0}, {"mean_incub", 5.42}, {"incub_period", 5.42}},
                           0., false, hsp::na::DOUBLE, false) {}

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
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Register Distributions  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

REGISTER_HOUSEPI_FUNCTION(Uniform_, ContinuosDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Normal_, ContinuosDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Lognormal_, ContinuosDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Beta_, ContinuosDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Gamma_, ContinuosDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(GammaMeanVar_, ContinuosDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Rayleigh_, ContinuosDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Weibull_, ContinuosDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Exponential_, ContinuosDistrRegistry)

REGISTER_HOUSEPI_FUNCTION(Poisson_, DiscreteDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Binomial_, DiscreteDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(NegativeBinomial_, DiscreteDistrRegistry)
REGISTER_HOUSEPI_FUNCTION(Geometric_, DiscreteDistrRegistry)

REGISTER_HOUSEPI_FUNCTION(InfecProfile_Ferretti20, InfecProfileRegistry)
