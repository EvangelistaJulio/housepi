/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS_MATH_DISTRIBUTIONS__H
#define HOUSEPI__TOOLS_MATH_DISTRIBUTIONS__H

#include <cmath>

// Suppress Boost warnings
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-copy"
#endif

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/geometric.hpp>

// Restore warnings
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#ifdef __clang__
#pragma clang diagnostic pop
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /// 1/sqrt(2pi)
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362105 /// 1/sqrt(2)
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /// log(sqrt(2pi))
#endif

// namespace containing numerical recipes used to compute distributions
namespace numr
{
    /// Logarithm of binomial coefficient.
    inline double log_binom_coeff(unsigned N, unsigned k)
    {
        if (N == 0 || k == N)
            return 0.;
        return std::lgamma(N + 1) - std::lgamma(k + 1) - std::lgamma(N - k + 1);
    }

    /// Log-ratio of binomial coefficients.
    inline double logratio_binom_coeff(unsigned N, unsigned ko, unsigned kn)
    {
        if (N == 0 || ko == kn)
            return 0;
        return std::lgamma(ko + 1) + std::lgamma(N - ko + 1) - std::lgamma(kn + 1) - std::lgamma(N - kn + 1);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace normal
    {
        inline double pdf(double x, double mean, double sd)
        {
            const double z = (x - mean) / sd;
            return M_1_SQRT_2PI * std::exp(-.5 * z * z) / sd;  
        }
        inline double cdf(double x, double mean, double sd)
        {
            return .5 * std::erfc(M_SQRT1_2 * (mean - x) / sd);
        }
        inline double logpdf(double x, double mean, double sd)
        {
            const double z = (x - mean) / sd;
            return -.5 * z * z - M_LN_SQRT_2PI - std::log(sd);
        }
        inline double logratio(double xo, double xn, double mean, double sd)
        {
            return (.5 * xo * xo + (xn - xo) * mean - .5 * xn * xn) / (sd * sd);
        }
        inline double quantile(double p, double mean, double sd)
        {
            static const boost::math::normal stdnorm(0., 1.);
            return mean + sd * boost::math::quantile(stdnorm, p);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  LOGNORMAL DISTRIBUTION METHODS  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace lognormal
    {
        inline double pdf(double x, double meanlog, double sdlog)
        {
            const double z = (std::log(x) - meanlog) / sdlog;
            return std::exp(-.5 * z * z) * M_1_SQRT_2PI / (x * sdlog);
        }
        inline double cdf(double x, double meanlog, double sdlog)
        {
            return .5 * std::erfc(M_SQRT1_2 * (meanlog - std::log(x)) / sdlog);
        }
        inline double logpdf(double x, double meanlog, double sdlog)
        {
            const double z = (std::log(x) - meanlog) / sdlog;
            return -.5 * z * z - std::log(x * sdlog) - M_LN_SQRT_2PI;
        }
        inline double logratio(double xo, double xn, double meanlog, double sdlog)
        {
            const double lxo = std::log(xo);
            const double lxn = std::log(xn);
            return .5 * ((lxo - meanlog) * lxo - (lxn - meanlog) * lxn) / (sdlog * sdlog) + lxo - lxn;
        }
        inline double quantile(double p, double meanlog, double sdlog)
        {
            static const boost::math::normal stdnorm(0., 1.);
            return std::exp(meanlog + sdlog * boost::math::quantile(stdnorm, p));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  BETA DISTRIBUTION METHODS  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace beta
    {
        inline double lbeta(double a, double b) { return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b); }

        inline double pdf(double x, double alpha, double beta)
        {
            return std::exp((alpha - 1.) * std::log(x) + (beta - 1.) * std::log1p(-x) - lbeta(alpha, beta));
        }
        inline double cdf(double x, double alpha, double beta)
        {
            return boost::math::ibeta(alpha, beta, x);
        }
        inline double logpdf(double x, double alpha, double beta)
        {
            return (alpha - 1.) * std::log(x) + (beta - 1.) * std::log1p(-x) - lbeta(alpha, beta);
        }
        inline double logratio(double xo, double xn, double alpha, double beta)
        {
            return (alpha - 1.) * std::log(xn / xo) + (beta - 1.) * std::log((1. - xn) / (1. - xo));
        }
        inline double quantile(double p, double alpha, double beta)
        {
            return boost::math::quantile(boost::math::beta_distribution<>(alpha, beta), p);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  GAMMA DISTRIBUTION METHODS  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace gamma
    {
        inline double pdf(double x, double shape, double scale)
        {
            return std::pow(x, shape - 1.) * std::exp(-shape * std::log(scale) - std::lgamma(shape) - x / scale);
        }
        inline double cdf(double x, double shape, double scale)
        {
            return boost::math::gamma_p(shape, x / scale);
        }
        inline double logpdf(double x, double shape, double scale)
        {
            const double t = x / scale;
            return shape * std::log(t) - std::log(x) - std::lgamma(shape) - t;
        }
        inline double logratio(double xo, double xn, double shape, double scale)
        {
            return (shape - 1.) * std::log(xn / xo) + (xo - xn) / scale;
        }
        inline double quantile(double p, double shape, double scale)
        {
            return boost::math::quantile(boost::math::gamma_distribution<>(shape, scale), p);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  RAYLEIGH DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace rayleigh
    {
        inline double pdf(double x, double scale)
        {
            const double scale2 = scale * scale;
            return x * std::exp(-.5 * x * x / scale2) / scale2;
        }
        inline double cdf(double x, double scale)
        {
            return 1. - std::exp(-.5 * x * x / (scale * scale));
        }
        inline double logpdf(double x, double scale)
        {
            const double scale2 = scale * scale;
            return std::log(x / scale2) - .5 * x * x / scale2;
        }
        inline double logratio(double xo, double xn, double scale)
        {
            return std::log(xn / xo) + .5 * (xo * xo - xn * xn) / (scale * scale);
        }
        inline double quantile(double p, double sigma)
        {
            return sigma * std::sqrt(-2. * std::log1p(-p));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  WEIBULL DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace weibull
    {
        inline double pdf(double x, double shape, double scale)
        {
            const double z = x / scale;
            return (shape / scale) * std::pow(z, shape - 1.) * std::exp(-std::pow(z, shape));
        }
        inline double cdf(double x, double shape, double scale)
        {
            return 1. - std::exp(-std::pow(x / scale, shape));
        }
        inline double logpdf(double x, double shape, double scale)
        {
            const double z = x / scale;
            return std::log(shape / scale) + (shape - 1.) * std::log(z) - std::pow(z, shape);
        }
        inline double logratio(double xo, double xn, double shape, double scale)
        {
            return (shape - 1.) * std::log(xn / xo) + std::pow(xo / scale, shape) - std::pow(xn / scale, shape);
        }
        inline double quantile(double p, double shape, double scale)
        {
            return scale * std::pow(-std::log1p(-p), 1. / shape);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  EXPONENTIAL DISTRIBUTION METHODS  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace exponential
    {
        inline double pdf(double x, double rate) { return rate * std::exp(-rate * x); }
        inline double cdf(double x, double rate) { return 1. - std::exp(-rate * x); }
        inline double logpdf(double x, double rate) { return std::log(rate) - rate * x; }
        inline double logratio(double xo, double xn, double rate) { return rate * (xo - xn); }
        inline double quantile(double p, double rate) { return -std::log1p(-p) / rate; }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  POISSON DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace poisson
    {
        inline double pdf(double x, double lambda)
        {
            return std::exp(x * std::log(lambda) - lambda - std::lgamma(x + 1.));
        }
        inline double cdf(double x, double lambda)
        {
            return 1. - boost::math::gamma_p(x + 1., lambda);
        }
        inline double logpdf(double x, double lambda)
        {
            return x * std::log(lambda) - lambda - std::lgamma(x + 1.);
        }
        inline double logratio(double xo, double xn, double lambda)
        {
            return (xn - xo) * std::log(lambda) + std::lgamma(xo + 1.) - std::lgamma(xn + 1.);
        }
        inline double quantile(double p, double lambda)
        {
            return boost::math::quantile(boost::math::poisson_distribution<>(lambda), p);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  BINOMIAL DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace binomial
    {
        inline double pdf(double x, double size, double prob)
        {
            return std::exp(log_binom_coeff(size, x) + x * std::log(prob) + (size - x) * std::log1p(-prob));
        }
        inline double cdf(double x, unsigned size, double prob)
        {
            if (x > size)
                return 1.;
            return boost::math::ibeta(size - x, x + 1, 1. - prob);
        }
        inline double logpdf(double x, double size, double prob)
        {
            return log_binom_coeff(size, x) + x * std::log(prob) + (size - x) * std::log1p(-prob);
        }
        inline double logratio(double xo, double xn, double size, double prob)
        {
            return logratio_binom_coeff(size, xo, xn) + (xn - xo) * std::log(prob / (1. - prob));
        }
        inline double quantile(double p, double size, double prob)
        {
            return boost::math::quantile(boost::math::binomial_distribution<>(size, prob), p);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  NEGATIVE BINOMIAL DISTRIBUTION METHODS  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace negbinomial
    {
        inline double pdf(double x, double size, double prob)
        {
            return std::exp(log_binom_coeff(x + size - 1., x) + x * std::log1p(-prob) + size * std::log(prob));
        }
        inline double cdf(double x, double size, double prob)
        {
            return boost::math::ibeta(size, x + 1., prob);
        }
        inline double logpdf(double x, double size, double prob)
        {
            return log_binom_coeff(x + size - 1., x) + x * std::log1p(-prob) + size * std::log(prob);
        }
        inline double logratio(double xo, double xn, double size, double prob)
        {
            return logratio_binom_coeff(size, xo, xn) + (xn - xo) * std::log1p(-prob);
        }
        inline double quantile(double p, double size, double prob)
        {
            return boost::math::quantile(boost::math::negative_binomial_distribution<>(size, prob), p);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  GEOMETRIC DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace geometric
    {
        inline double pdf(double x, double prob) { return prob * std::pow(1. - prob, x); }
        inline double cdf(double x, double prob) { return 1. - std::pow(1. - prob, x + 1); }
        inline double logpdf(double x, double prob) { return std::log(prob) + x * std::log1p(-prob); }
        inline double logratio(double xo, double xn, double prob) { return (xn - xo) * std::log1p(-prob); }
        inline double quantile(double p, double prob)
        {
            return boost::math::quantile(boost::math::geometric_distribution<>(prob), p);
        }
    }
} // NUMR NAMESPACE

#endif // HOUSEPI__TOOLS_MATH_DISTRIBUTIONS__H
