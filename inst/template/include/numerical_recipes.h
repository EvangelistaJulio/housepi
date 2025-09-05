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

#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <string>

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /// 1/sqrt(2pi)
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362105 /// 1/sqrt(2)
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /// log(sqrt(2pi))
#endif

#ifndef HOUSEPI__NUMERICAL_RECIPES__H
#define HOUSEPI__NUMERICAL_RECIPES__H

// namespace containing numerical recipes used to compute distributions
namespace numr
{
    constexpr unsigned MAX_ITER = 1000;
    constexpr double EPSILON = 1.e-12;
    constexpr double MIN_DBL = 1.e-30;

    using std::exp;
    using std::lgamma;
    using std::log;
    using std::pow;

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
    unsigned long binom_coeff(unsigned short N, unsigned short k);

    /// Logarithm of binomial coefficient.
    double log_binom_coeff(unsigned short N, unsigned short k);

    /// Log-ratio of binomial coefficients.
    double logratio_binom_coeff(unsigned short N, unsigned short ko, unsigned short kn);

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  REGULARISED INCOMPLETE GAMMA FUNCTION  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// UPDATE: To better reference this:
    /// Press, Teukolsky, Vetterling, and Flannery
    /// Numerical Recipes 3rd Edition: The Art of Scientific Computing
    namespace gamma_p_impl
    {
        /// Lower regulairzed incomplete gamma function using series expansion.
        double series(double shape, double z);

        /// Upper regularized incomplete gamma function using continued fraction.
        double continued_fraction(double shape, double z);
    }

    /// Regularized incomplete gamma function.
    double gamma_p(double shape, double z);

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  REGULARISED INCOMPLETE BETA FUNCTION  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Logarithm of the beta function.
    inline double lbeta(double a, double b)
    {
        return lgamma(a) + lgamma(b) - lgamma(a + b);
    }

    /// UPDATE: to better reference this
    /// (C) Copr. 1986-92 Numerical Recipes Software Vn2-5.
    namespace ibeta_impl
    {
        double continued_fraction(double a, double b, double x);
    }

    /// Regularized incomplete beta function.
    double ibeta(double a, double b, double x);

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace normal
    {
        inline double pdf(double x, double mean, double sd)
        {
            const double z = (x - mean) / sd;
            return M_1_SQRT_2PI * exp(-.5 * z * z) / sd;
        }
        inline double cdf(double x, double mean, double sd)
        {
            return .5 + .5 * std::erf(M_SQRT1_2 * (x - mean) / sd);
        }
        inline double logpdf(double x, double mean, double sd)
        {
            const double z = (x - mean) / sd;
            return -.5 * z * z - M_LN_SQRT_2PI - log(sd);
        }
        inline double logratio(double xo, double xn, double mean, double sd)
        {
            return (.5 * xo * xo + (xn - xo) * mean - .5 * xn * xn) / (sd * sd);
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
            return exp(-.5 * z * z) * M_1_SQRT_2PI / (x * sdlog);
        }
        inline double cdf(double x, double meanlog, double sdlog)
        {
            return .5 + .5 * std::erf(M_SQRT1_2 * (log(x) - meanlog) / sdlog);
        }
        inline double logpdf(double x, double meanlog, double sdlog)
        {
            const double z = (log(x) - meanlog) / sdlog;
            return -.5 * z * z - log(x * sdlog) - M_LN_SQRT_2PI;
        }
        inline double logratio(double xo, double xn, double meanlog, double sdlog)
        {
            const double lxo = log(xo);
            const double lxn = log(xn);
            return .5 * ((lxo - meanlog) * lxo - (lxn - meanlog) * lxn) / (sdlog * sdlog) + lxo - lxn;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  BETA DISTRIBUTION METHODS  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace beta
    {
        inline double pdf(double x, double shape1, double shape2)
        {
            return exp((shape1 - 1.) * log(x) + (shape2 - 1.) * log(1. - x) - lbeta(shape1, shape2));
        }
        inline double cdf(double x, double shape1, double shape2)
        {
            return ibeta(shape1, shape2, x);
        }
        inline double logpdf(double x, double shape1, double shape2)
        {
            return (shape1 - 1.) * log(x) + (shape2 - 1.) * log(1. - x) - lbeta(shape1, shape2);
        }
        inline double logratio(double xo, double xn, double shape1, double shape2)
        {
            return (shape1 - 1.) * log(xn / xo) + (shape2 - 1.) * log((1. - xn) / (1. - xo));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  GAMMA DISTRIBUTION METHODS  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace gamma
    {
        inline double pdf(double x, double shape, double scale)
        {
            return pow(x, shape - 1) * exp(-shape * log(scale) - lgamma(shape) - x / scale);
        }
        inline double cdf(double x, double shape, double scale)
        {
            return gamma_p(shape, x / scale);
        }
        inline double logpdf(double x, double shape, double scale)
        {
            const double t = x / scale;
            return shape * log(t) - log(x) - lgamma(shape) - t;
        }
        inline double logratio(double xo, double xn, double shape, double scale)
        {
            return (shape - 1) * log(xn / xo) + (xo - xn) / scale;
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
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  EXPONENTIAL DISTRIBUTION METHODS  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace exponential
    {
        /// @param x
        /// @param r rate
        inline double pdf(double x, double r)
        {
            return r * exp(-r * x);
        }
        inline double cdf(double x, double r)
        {
            return 1. - exp(-r * x);
        }
        inline double logpdf(double x, double r)
        {
            return log(r) - r * x;
        }
        inline double logratio(double xo, double xn, double r) { return r * (xo - xn); }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  POISSON DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace poisson
    {
        /// @param x
        /// @param L lambda
        inline double pdf(double x, double L)
        {
            return exp(x * log(L) - L - lgamma(x + 1.));
        }
        inline double cdf(double x, double L)
        {
            return 1. - gamma_p(x + 1., L);
        }
        inline double logpdf(double x, double L)
        {
            return x * log(L) - L - lgamma(x + 1.);
        }
        inline double logratio(double xo, double xn, double L)
        {
            return (xn - xo) * log(L) + lgamma(xo + 1.) - lgamma(xn + 1.);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  BINOMIAL DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace binomial
    {
        /// @param x
        /// @param s size
        /// @param p probability
        inline double pdf(double x, double s, double p)
        {
            return exp(log_binom_coeff(s, x) + x * log(p) + (s - x) * log(1. - p));
        }
        inline double cdf(double x, double s, double p)
        {
            return ibeta(s - x, x + 1, 1. - p);
        }
        inline double logpdf(double x, double s, double p)
        {
            return log_binom_coeff(s, x) + x * log(p) + (s - x) * log(1. - p);
        }
        inline double logratio(double xo, double xn, double s, double p)
        {
            return logratio_binom_coeff(s, xo, xn) + (xn - xo) * log(p / (1. - p));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  NEGATIVE BINOMIAL DISTRIBUTION METHODS  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace negbinomial
    {
        /// @param x
        /// @param s size
        /// @param p probability
        inline double pdf(double x, double s, double p)
        {
            return exp(log_binom_coeff(x + s - 1., x) + x * log(1. - p) + s * log(p));
        }
        inline double cdf(double x, double s, double p)
        {
            return ibeta(s, x + 1., p);
        }
        inline double logpdf(double x, double s, double p)
        {
            return log_binom_coeff(x + s - 1., x) + x * log(1. - p) + s * log(p);
        }
        inline double logratio(double xo, double xn, double s, double p)
        {
            double output = 0.;
            // logratio of the binomial coefficients
            for (double i = 1.; i < s; ++i)
                output += log((i + xn) / (i + xo));
            return output + (xn - xo) * log(1. - p);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  GEOMETRIC DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace geometric
    {
        /// @param x
        /// @param p probability
        inline double pdf(double x, double p)
        {
            return p * pow(1. - p, x);
        }
        inline double cdf(double x, double p)
        {
            return 1. - pow(1. - p, x + 1);
        }
        inline double logpdf(double x, double p)
        {
            return log(p) + x * log(1. - p);
        }
        inline double logratio(double xo, double xn, double p)
        {
            return (xn - xo) * log(1. - p);
        }
    }
} // NUMR NAMESPACE

#endif // HOUSEPI__NUMERICAL_RECIPES__H
