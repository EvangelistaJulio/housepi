/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * numerical_recipes.h
 *
 */

#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "tools.h"

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362105 /* 1/sqrt(2) */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /* log(sqrt(2pi)) */
#endif

#ifndef HOUSEPI__NUMERICAL_RECIPES__H
#define HOUSEPI__NUMERICAL_RECIPES__H

// namespace containing numerical recipes used to compute distributions
namespace numr
{
    constexpr unsigned MAX_ITER = 1000;
    constexpr double EPSILON = 1.e-12;
    constexpr double MIN_DBL = 1.e-30;

    using ushort = unsigned short int;
    using ulong = unsigned long int;

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  COMBINATORICS FUNCTIONS  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* Standard definition for the factorial. */
    inline ulong factorial(ulong n)
    {
        ulong result = 1;
        for (ulong i = 2; i <= n; ++i)
            result *= i;
        return result;
    }

    /* Binomial coefficient of N elements in k positions */
    inline ulong binom_coeff(ushort N, ushort k)
    {
        if ((N > 0) && (k != N))
        {
            ulong result = N;
            ushort stop = (k > N - k) ? k : N - k; // use symmetry
            for (ushort i = 1; i < stop; ++i)
                result = (result * (N - i)) / (i + 1);
            return result;
        }

        return 1;
    }

    /* Logarithm of binomial coefficient. */
    inline double log_binom_coeff(ushort N, ushort k)
    {
        if ((N > 0) && (k != N))
        {
            double result = std::log(N);
            ushort stop = (k > N - k) ? k : N - k; // use symmetry
            for (ushort i = 1; i < stop; ++i)
                result += std::log(N - i) - std::log(i + 1);
            return result;
        }
        return 0.;
    }

    /* Log-ratio of binomial coefficients. */
    inline double logratio_binom_coeff(ushort N, ushort kold, ushort knew)
    {
        if ((N > 0) && (kold != knew))
        {
            double result = 0.;
            if (kold > knew)
                std::swap(kold, knew);

            for (ushort i = kold + 1; i <= knew; ++i)
                result += std::log(N - i + 1) - std::log(i);
            return result;
        }
        return 0.;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  REGULARISED INCOMPLETE GAMMA FUNCTION  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /*
     * UPDATE: To better reference this:
     * Press, Teukolsky, Vetterling, and Flannery
     * Numerical Recipes 3rd Edition: The Art of Scientific Computing
     */
    namespace gamma_p_impl
    {
        /* Lower regulairzed incomplete gamma function using series expansion. */
        inline double series(double shape, double z)
        {
            double sum = 1. / shape;
            double term = sum;
            for (unsigned m = 1; m <= MAX_ITER; ++m)
            {
                term *= z / (shape + m);
                sum += term;

                if (term < EPSILON * sum)
                    return sum * std::exp(-z + shape * std::log(z) - std::lgamma(shape));
            }

            throw std::overflow_error(
                "Overflow in gamma_p function for x/scale=" + std::to_string(z) +
                ", shape=" + std::to_string(shape) + ".");
        }

        /* Upper regularized incomplete gamma function using continued fraction. */
        inline double continued_fraction(double shape, double z)
        {
            double b = z + 1. - shape;
            double c = 1. / MIN_DBL;
            double d = 1. / b;
            double h = d;
            double delta = 0.;

            for (unsigned m = 1; m <= MAX_ITER; ++m)
            {
                double an = (shape - m) * m;
                b += 2.;
                d = an * d + b;
                if (std::fabs(d) < MIN_DBL)
                    d = MIN_DBL;
                c = b + an / c;
                if (std::fabs(c) < MIN_DBL)
                    c = MIN_DBL;
                d = 1. / d;
                delta = c * d;
                h *= delta;

                if (std::fabs(delta - 1.) < EPSILON)
                    return 1. - std::exp(-z + shape * std::log(z) - std::lgamma(shape)) * h;
            }

            throw std::overflow_error(
                "Overflow in gamma_p function for x/scale=" + std::to_string(z) +
                ", shape=" + std::to_string(shape) + ".");
        }
    }

    /* Regularized incomplete gamma function. */
    inline double gamma_p(double shape, double z)
    {
        if (z < shape + 10.)
            return gamma_p_impl::series(shape, z);
        else
            return gamma_p_impl::continued_fraction(shape, z);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////  REGULARISED INCOMPLETE BETA FUNCTION  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /* (C) Copr. 1986-92 Numerical Recipes Software Vn2-5. */
    namespace ibeta_impl
    {
        inline double continued_fraction(double a, double b, double x)
        {
            const double qab = a + b,
                         qap = a + 1.,
                         qam = a - 1.;
            double c = 1.,
                   d = 1. - qab * x / qap;

            if (std::fabs(d) < MIN_DBL)
                d = MIN_DBL;
            d = 1. / d;

            double h = d;
            for (unsigned m = 1; m <= MAX_ITER; ++m)
            {
                unsigned m2 = 2 * m;
                double aa = (b - m) * m * x / ((qam + m2) * (a + m2));

                d = 1. + aa * d;
                if (std::fabs(d) < MIN_DBL)
                    d = MIN_DBL;

                c = 1. + aa / c;
                if (std::fabs(c) < MIN_DBL)
                    c = MIN_DBL;

                d = 1. / d;
                h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));

                d = 1. + aa * d;
                if (std::fabs(d) < MIN_DBL)
                    d = MIN_DBL;

                c = 1. + aa / c;
                if (std::fabs(c) < MIN_DBL)
                    c = MIN_DBL;

                d = 1. / d;

                double del = d * c;
                h *= del;

                if (std::fabs(del - 1.) < EPSILON)
                    return h;
            }

            throw std::overflow_error(
                "Overflow in ibeta function for x=" + std::to_string(x) +
                ", a=" + std::to_string(a) + ", b=" + std::to_string(b) + ".");
        }
    }

    /* Logarithm of the beta function. */
    inline double lbeta(double a, double b)
    {
        return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
    }

    /* Regularized incomplete beta function. */
    inline double ibeta(double a, double b, double x)
    {
        if (x == 0.)
            return 0.;
        else if (x == 1.)
            return 1.;
        else
        {
            const double lbt = a * std::log(x) + b * std::log(1.0 - x) - lbeta(a, b);

            if (x < (a + 1.) / (a + b + 2.))
                return std::exp(lbt + std::log(ibeta_impl::continued_fraction(a, b, x) / a));
            else
                return 1. - std::exp(lbt + std::log(ibeta_impl::continued_fraction(b, a, 1. - x) / b));
        }
    }

    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------------------//

    /*
     * Methods used for the definition of the derived distribution classes.
     *
     */

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace normal
    {
        inline double pdf(double x, double mean, double sd)
        {
            const double ptm = (x - mean) / sd;
            return M_1_SQRT_2PI * std::exp(-.5 * ptm * ptm) / sd;
        }
        inline double cdf(double x, double mean, double sd)
        {
            return .5 + .5 * std::erf(M_SQRT1_2 * (x - mean) / sd);
        }
        inline double logpdf(double x, double mean, double sd)
        {
            const double ptm = (x - mean) / sd;
            return -.5 * ptm * ptm - M_LN_SQRT_2PI - std::log(sd);
        }
        inline double logratio(double xold, double xnew, double mean, double sd)
        {
            return .5 * (xold * xold + (xnew - xold) * mean - xnew * xnew) / (sd * sd);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  LOGNORMAL DISTRIBUTION METHODS  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace lognormal
    {
        inline double pdf(double x, double meanlog, double sdlog)
        {
            const double ptm = (std::log(x) - meanlog) / sdlog;
            return std::exp(-.5 * ptm * ptm) * M_1_SQRT_2PI / (x * sdlog);
        }
        inline double cdf(double x, double meanlog, double sdlog)
        {
            return .5 + .5 * std::erf(M_SQRT1_2 * (std::log(x) - meanlog) / sdlog);
        }
        inline double logpdf(double x, double meanlog, double sdlog)
        {
            const double ptm = (std::log(x) - meanlog) / sdlog;
            return -.5 * ptm * ptm - std::log(x * sdlog) - M_LN_SQRT_2PI;
        }
        inline double logratio(double xold, double xnew, double meanlog, double sdlog)
        {
            const double logxold = std::log(xold);
            const double logxnew = std::log(xnew);
            return .5 * ((logxold - meanlog) * logxold - (logxnew - meanlog) * logxnew) / (sdlog * sdlog) + logxold - logxnew;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  BETA DISTRIBUTION METHODS  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace beta
    {
        inline double pdf(double x, double shape1, double shape2)
        {
            return std::exp((shape1 - 1.) * std::log(x) + (shape2 - 1.) * std::log(1. - x) - lbeta(shape1, shape2));
        }
        inline double cdf(double x, double shape1, double shape2)
        {
            return ibeta(shape1, shape2, x);
        }
        inline double logpdf(double x, double shape1, double shape2)
        {
            return (shape1 - 1.) * std::log(x) + (shape2 - 1.) * std::log(1. - x) - lbeta(shape1, shape2);
        }
        inline double logratio(double xold, double xnew, double shape1, double shape2)
        {
            return (shape1 - 1.) * std::log(xnew / xold) + (shape2 - 1.) * std::log((1. - xnew) / (1. - xold));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  GAMMA DISTRIBUTION METHODS  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace gamma
    {
        inline double pdf(double x, double shape, double scale)
        {
            return std::pow(x, shape - 1) * std::exp(-shape * std::log(scale) - std::lgamma(shape) - x / scale);
        }
        inline double cdf(double x, double shape, double scale)
        {
            return gamma_p(shape, x / scale);
        }
        inline double logpdf(double x, double shape, double scale)
        {
            const double t = x / scale;
            return shape * std::log(t) - std::log(x) - std::lgamma(shape) - t;
        }
        inline double logratio(double xold, double xnew, double shape, double scale)
        {
            return (shape - 1) * std::log(xnew / xold) + (xold - xnew) / scale;
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
        inline double logratio(double xold, double xnew, double scale)
        {
            return std::log(xnew / xold) + .5 * (xold * xold - xnew * xnew) / (scale * scale);
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
        inline double logratio(double xold, double xnew, double shape, double scale)
        {
            return (shape - 1.) * std::log(xnew / xold) + std::pow(xold / scale, shape) - std::pow(xnew / scale, shape);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  EXPONENTIAL DISTRIBUTION METHODS  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace exponential
    {
        inline double pdf(double x, double rate)
        {
            return rate * std::exp(-rate * x);
        }
        inline double cdf(double x, double rate)
        {
            return 1. - std::exp(-rate * x);
        }
        inline double logpdf(double x, double rate)
        {
            return std::log(rate) - rate * x;
        }
        inline double logratio(double xold, double xnew, double rate)
        {
            return rate * (xold - xnew);
        }
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
            return 1. - gamma_p(x + 1., lambda);
        }
        inline double logpdf(double x, double lambda)
        {
            return x * std::log(lambda) - lambda - std::lgamma(x + 1.);
        }
        inline double logratio(double xold, double xnew, double lambda)
        {
            return (xnew - xold) * std::log(lambda) + std::lgamma(xold + 1.) - std::lgamma(xnew + 1.);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  BINOMIAL DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace binomial
    {
        inline double pdf(double x, double size, double prob)
        {
            return std::exp(log_binom_coeff(size, x) + x * std::log(prob) + (size - x) * std::log(1. - prob));
        }
        inline double cdf(double x, double size, double prob)
        {
            return ibeta(size - x, x + 1, 1. - prob);
        }
        inline double logpdf(double x, double size, double prob)
        {
            return log_binom_coeff(size, x) + x * std::log(prob) + (size - x) * std::log(1. - prob);
        }
        inline double logratio(double xold, double xnew, double size, double prob)
        {
            return logratio_binom_coeff(size, xold, xnew) + (xnew - xold) * std::log(prob / (1. - prob));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  NEGATIVE BINOMIAL DISTRIBUTION METHODS  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace negbinomial
    {
        inline double pdf(double x, double size, double prob)
        {
            return std::exp(log_binom_coeff(x + size - 1, x) + x * std::log(1. - prob) + size * std::log(prob));
        }
        inline double cdf(double x, double size, double prob)
        {
            return ibeta(size, x + 1, prob);
        }
        inline double logpdf(double x, double size, double prob)
        {
            return log_binom_coeff(x + size - 1, x) + x * std::log(1. - prob) + size * std::log(prob);
        }
        inline double logratio(double xold, double xnew, double size, double prob)
        {
            double output = 0.;
            // logratio of the binomial coefficients
            for (double i = 1.; i < size; ++i)
                output += std::log((i + xnew) / (i + xold));
            return output + (xnew - xold) * std::log(1. - prob);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////  GEOMETRIC DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace geometric
    {
        inline double pdf(double x, double proba)
        {
            return proba * std::pow(1. - proba, x);
        }
        inline double cdf(double x, double proba)
        {
            return 1. - std::pow(1. - proba, x + 1);
        }
        inline double logpdf(double x, double proba)
        {
            return std::log(proba) + x * std::log(1. - proba);
        }
        inline double logratio(double xold, double xnew, double proba)
        {
            return (xnew - xold) * std::log(1. - proba);
        }
    }
} // NUMR NAMESPACE

#endif // HOUSEPI__NUMERICAL_RECIPES__H
