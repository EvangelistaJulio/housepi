/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * numerical_recipes.h
 *
 */

#ifndef HOUSEPI__NUMERICAL_RECIPES__H
#define HOUSEPI__NUMERICAL_RECIPES__H

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
    using std::exp;
    using std::lgamma;
    using std::log;
    using std::log1p;
    using std::pow;
    using std::sqrt;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  COMBINATORICS FUNCTIONS  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Standard definition for the factorial.
    inline unsigned long factorial(unsigned long n)
    {
        unsigned long result = 1;
        for (unsigned long i = 2; i <= n; ++i)
            result *= i;
        return result;
    }

    /// Binomial coefficient of N elements in k positions
    inline unsigned long binom_coeff(unsigned short N, unsigned short k)
    {
        if (N == 0 || k == N)
            return 1;
        // note: use symmetry
        if (k > N - k)
            k = N - k;
        unsigned long result = 1;
        for (unsigned short i = 1; i <= k; ++i)
            result = result * (N - k + i) / i;
        return result;
    }

    /// Logarithm of binomial coefficient.
    inline double log_binom_coeff(unsigned short N, unsigned short k)
    {
        if (N == 0 || k == N)
            return 0.;
        return lgamma(N + 1) - lgamma(k + 1) - lgamma(N - k + 1);
    }

    /// Log-ratio of binomial coefficients.
    inline double logratio_binom_coeff(unsigned short N, unsigned short ko, unsigned short kn)
    {
        if (N == 0 || ko == kn)
            return 0;
        return lgamma(ko + 1) + lgamma(N - ko + 1) - lgamma(kn + 1) - lgamma(N - kn + 1);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace normal
    {
        inline double pdf(double x, double mean, double sd)
        {
            const double z = (x - mean) / sd;
            return static_cast<double>(M_1_SQRT_2PI) * exp(-.5 * z * z) / sd;
        }
        inline double cdf(double x, double mean, double sd)
        {
            return .5 + .5 * std::erf(M_SQRT1_2 * (x - mean) / sd);
        }
        inline double logpdf(double x, double mean, double sd)
        {
            const double z = (x - mean) / sd;
            return -.5 * z * z - static_cast<double>(M_LN_SQRT_2PI) - log(sd);
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
            const double z = (log(x) - meanlog) / sdlog;
            return exp(-.5 * z * z) * static_cast<double>(M_1_SQRT_2PI) / (x * sdlog);
        }
        inline double cdf(double x, double meanlog, double sdlog)
        {
            return .5 + .5 * std::erf(static_cast<double>(M_SQRT1_2) * (log(x) - meanlog) / sdlog);
        }
        inline double logpdf(double x, double meanlog, double sdlog)
        {
            const double z = (log(x) - meanlog) / sdlog;
            return -.5 * z * z - log(x * sdlog) - static_cast<double>(M_LN_SQRT_2PI);
        }
        inline double logratio(double xo, double xn, double meanlog, double sdlog)
        {
            const double lxo = log(xo);
            const double lxn = log(xn);
            return .5 * ((lxo - meanlog) * lxo - (lxn - meanlog) * lxn) / (sdlog * sdlog) + lxo - lxn;
        }
        inline double quantile(double p, double meanlog, double sdlog)
        {
            static const boost::math::normal stdnorm(0., 1.);
            return exp(meanlog + sdlog * boost::math::quantile(stdnorm, p));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  BETA DISTRIBUTION METHODS  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace beta
    {
        inline double lbeta(double a, double b) { return lgamma(a) + lgamma(b) - lgamma(a + b); }

        inline double pdf(double x, double alpha, double beta)
        {
            return exp((alpha - 1.) * log(x) + (beta - 1.) * log1p(-x) - lbeta(alpha, beta));
        }
        inline double cdf(double x, double alpha, double beta)
        {
            return boost::math::ibeta(alpha, beta, x);
        }
        inline double logpdf(double x, double alpha, double beta)
        {
            return (alpha - 1.) * log(x) + (beta - 1.) * log1p(-x) - lbeta(alpha, beta);
        }
        inline double logratio(double xo, double xn, double alpha, double beta)
        {
            return (alpha - 1.) * log(xn / xo) + (beta - 1.) * log((1. - xn) / (1. - xo));
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
            return pow(x, shape - 1.) * exp(-shape * log(scale) - lgamma(shape) - x / scale);
        }
        inline double cdf(double x, double shape, double scale)
        {
            return boost::math::gamma_p(shape, x / scale);
        }
        inline double logpdf(double x, double shape, double scale)
        {
            const double t = x / scale;
            return shape * log(t) - log(x) - lgamma(shape) - t;
        }
        inline double logratio(double xo, double xn, double shape, double scale)
        {
            return (shape - 1.) * log(xn / xo) + (xo - xn) / scale;
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
            return x * exp(-.5 * x * x / scale2) / scale2;
        }
        inline double cdf(double x, double scale)
        {
            return 1. - exp(-.5 * x * x / (scale * scale));
        }
        inline double logpdf(double x, double scale)
        {
            const double scale2 = scale * scale;
            return log(x / scale2) - .5 * x * x / scale2;
        }
        inline double logratio(double xo, double xn, double scale)
        {
            return log(xn / xo) + .5 * (xo * xo - xn * xn) / (scale * scale);
        }
        inline double quantile(double p, double sigma)
        {
            return sigma * sqrt(-2. * log1p(-p));
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
            return (shape / scale) * pow(z, shape - 1.) * exp(-pow(z, shape));
        }
        inline double cdf(double x, double shape, double scale)
        {
            return 1. - exp(-pow(x / scale, shape));
        }
        inline double logpdf(double x, double shape, double scale)
        {
            const double z = x / scale;
            return log(shape / scale) + (shape - 1.) * log(z) - pow(z, shape);
        }
        inline double logratio(double xo, double xn, double shape, double scale)
        {
            return (shape - 1.) * log(xn / xo) + pow(xo / scale, shape) - pow(xn / scale, shape);
        }
        inline double quantile(double p, double shape, double scale)
        {
            return scale * pow(-log1p(-p), 1. / shape);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  EXPONENTIAL DISTRIBUTION METHODS  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace exponential
    {
        inline double pdf(double x, double rate) { return rate * exp(-rate * x); }
        inline double cdf(double x, double rate) { return 1. - exp(-rate * x); }
        inline double logpdf(double x, double rate) { return log(rate) - rate * x; }
        inline double logratio(double xo, double xn, double rate) { return rate * (xo - xn); }
        inline double quantile(double p, double rate) { return -log1p(-p) / rate; }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  POISSON DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace poisson
    {
        inline double pdf(double x, double lambda)
        {
            return exp(x * log(lambda) - lambda - lgamma(x + 1.));
        }
        inline double cdf(double x, double lambda)
        {
            return 1. - boost::math::gamma_p(x + 1., lambda);
        }
        inline double logpdf(double x, double lambda)
        {
            return x * log(lambda) - lambda - lgamma(x + 1.);
        }
        inline double logratio(double xo, double xn, double lambda)
        {
            return (xn - xo) * log(lambda) + lgamma(xo + 1.) - lgamma(xn + 1.);
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
            return exp(log_binom_coeff(size, x) + x * log(prob) + (size - x) * log1p(-prob));
        }
        inline double cdf(double x, unsigned size, double prob)
        {
            if (x > size)
                return 1.;
            return boost::math::ibeta(size - x, x + 1, 1. - prob);
        }
        inline double logpdf(double x, double size, double prob)
        {
            return log_binom_coeff(size, x) + x * log(prob) + (size - x) * log1p(-prob);
        }
        inline double logratio(double xo, double xn, double size, double prob)
        {
            return logratio_binom_coeff(size, xo, xn) + (xn - xo) * log(prob / (1. - prob));
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
            return exp(log_binom_coeff(x + size - 1., x) + x * log1p(-prob) + size * log(prob));
        }
        inline double cdf(double x, double size, double prob)
        {
            return boost::math::ibeta(size, x + 1., prob);
        }
        inline double logpdf(double x, double size, double prob)
        {
            return log_binom_coeff(x + size - 1., x) + x * log1p(-prob) + size * log(prob);
        }
        inline double logratio(double xo, double xn, double size, double prob)
        {
            double output = 0.;
            // logratio of the binomial coefficients
            for (double i = 1.; i < size; ++i)
                output += log((i + xn) / (i + xo));
            return output + (xn - xo) * log1p(-prob);
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
        inline double pdf(double x, double prob) { return prob * pow(1. - prob, x); }
        inline double cdf(double x, double prob) { return 1. - pow(1. - prob, x + 1); }
        inline double logpdf(double x, double prob) { return log(prob) + x * log1p(-prob); }
        inline double logratio(double xo, double xn, double prob) { return (xn - xo) * log1p(-prob); }
        inline double quantile(double p, double prob)
        {
            return boost::math::quantile(boost::math::geometric_distribution<>(prob), p);
        }
    }
} // NUMR NAMESPACE

#endif // HOUSEPI__NUMERICAL_RECIPES__H
