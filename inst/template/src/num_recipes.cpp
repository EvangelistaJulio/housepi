/*
 * num_recipes.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include <cmath>
#include <algorithm>

#include "num_recipes.h"

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /* 1/sqrt(2pi) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362105 /* 1/sqrt(2) */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /* log(sqrt(2pi)) */
#endif

// namespace for numerical recipes used for distributions
namespace numr
{
    constexpr int MAX_ITER = 1000;
    constexpr double EPSILON = 1.e-12;
    constexpr double MIN_DBL = 1.e-30;

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
    ulong factorial(ulong n)
    {
        ulong result = 1;
        for (ulong i = 2; i <= n; ++i)
            result *= i;
        return result;
    }

    /* Binomial coefficient of N elements in k positions */
    ulong binom_coeff(ushort N, ushort k)
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
    double log_binom_coeff(ushort N, ushort k)
    {
        if ((N > 0) && (k != N))
        {
            double result = log(N);
            ushort stop = (k > N - k) ? k : N - k; // use symmetry
            for (ushort i = 1; i < stop; ++i)
                result += log(N - i) - log(i + 1);
            return result;
        }
        return 0.;
    }

    /* Log-ratio of binomial coefficients. */
    double logratio_binom_coeff(ushort N, ushort kold, ushort knew)
    {
        if ((N > 0) && (kold != knew))
        {
            double result = 0.;
            if (kold > knew)
                std::swap(kold, knew);

            for (ushort i = kold + 1; i <= knew; ++i)
                result += log(N - i + 1) - log(i);
            return result;
        }
        return 0.;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  REGULARISED INCOMPLETE GAMMA FUNCTION  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    // Note: I could use Rmath.h but with limits in portability.
    /*
    #ifdef _HOUSEPI_RCPP
    #include <Rmath.h>
    // ibeta : regularized incomplete beta function.
    double ibeta(double a, double b, double x) { return R::pbeta(x, a, b, true, false); }

    // gamma_p : regularized incomplete beta function.
    double gamma_p(double shape, double z) { return R::pgamma(z, shape, 1., true, false); }
    #endif
    */

    /*
     * FIXME: To better reference this:
     * Press, Teukolsky, Vetterling, and Flannery
     * Numerical Recipes 3rd Edition: The Art of Scientific Computing
     */
    namespace gamma_p_impl
    {
        /* Lower regulairzed incomplete gamma function using series expansion. */
        double series(double shape, double z)
        {
            double sum = 1. / shape;
            double term = sum;
            for (int m = 1; m <= MAX_ITER; ++m)
            {
                term *= z / (shape + m);
                sum += term;

                if (term < EPSILON * sum)
                    return sum * std::exp(-z + shape * std::log(z) - std::lgamma(shape));
            }

            return std::numeric_limits<double>::quiet_NaN();
        }

        /* Upper regularized incomplete gamma function using continued fraction. */
        double continued_fraction(double shape, double z)
        {
            double b = z + 1. - shape;
            double c = 1. / MIN_DBL;
            double d = 1. / b;
            double h = d;
            double delta = 0.;

            for (int m = 1; m <= MAX_ITER; ++m)
            {
                double an = m * (shape - m);
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

            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    /* Regularized incomplete gamma function. */
    double gamma_p(double shape, double z)
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
        double continued_fraction(double a, double b, double x)
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
            for (int m = 1; m <= MAX_ITER; ++m)
            {
                int m2 = 2 * m;
                double aa = m * (b - m) * x / ((qam + m2) * (a + m2));

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

                if (fabs(del - 1.) < EPSILON)
                    return h;
            }

            return std::numeric_limits<double>::quiet_NaN();
        }
    }

    /* Logarithm of the beta function. */
    double lbeta(double a, double b)
    {
        return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
    }

    /* Regularized incomplete beta function. */
    double ibeta(double a, double b, double x)
    {
        if (x == 0.)
            return 0.;
        else if (x == 1.)
            return 1.;
        else
        {
            const double lbt = a * log(x) + b * log(1.0 - x) - lbeta(a, b);

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
        double pdf(double x, double mean, double sd)
        {
            const double ptm = (x - mean) / sd;
            return M_1_SQRT_2PI * std::exp(-.5 * ptm * ptm) / sd;
        }
        double cdf(double x, double mean, double sd)
        {
            return .5 + .5 * std::erf(M_SQRT1_2 * (x - mean) / sd);
        }
        double logpdf(double x, double mean, double sd)
        {
            const double ptm = (x - mean) / sd;
            return -.5 * ptm * ptm - M_LN_SQRT_2PI - std::log(sd);
        }
        double logratio(double xold, double xnew, double mean, double sd)
        {
            return .5 * (xold * xold + (xnew - xold) * mean - xnew * xnew) / (sd * sd);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  LOGNORMAL DISTRIBUTION METHODS  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace lognormal
    {
        double pdf(double x, double meanlog, double sdlog)
        {
            if (x <= 0.)
                return 0.;
            else
            {
                const double ptm = (std::log(x) - meanlog) / sdlog;
                return std::exp(-.5 * ptm * ptm) * M_1_SQRT_2PI / (x * sdlog);
            }
        }
        double cdf(double x, double meanlog, double sdlog)
        {
            if (x <= 0.)
                return 0.;
            else
                return .5 + .5 * std::erf(M_SQRT1_2 * (std::log(x) - meanlog) / sdlog);
        }
        double logpdf(double x, double meanlog, double sdlog)
        {
            if (x <= 0.)
                return std::numeric_limits<double>::quiet_NaN();
            else
            {
                const double ptm = (std::log(x) - meanlog) / sdlog;
                return -.5 * ptm * ptm - std::log(x * sdlog) - M_LN_SQRT_2PI;
            }
        }
        double logratio(double xold, double xnew, double meanlog, double sdlog)
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
        double pdf(double x, double shape1, double shape2)
        {
            return std::exp((shape1 - 1.) * std::log(x) + (shape2 - 1.) * std::log(1. - x) - lbeta(shape1, shape2));
        }
        double cdf(double x, double shape1, double shape2)
        {
            return ibeta(shape1, shape2, x);
        }
        double logpdf(double x, double shape1, double shape2)
        {
            return (shape1 - 1.) * std::log(x) + (shape2 - 1.) * std::log(1. - x) - lbeta(shape1, shape2);
        }
        double logratio(double xold, double xnew, double shape1, double shape2)
        {
            return (shape1 - 1.) * std::log(xnew / xold) + (shape2 - 1.) * std::log((1. - xnew) / (1. - xold));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  GAMMA DISTRIBUTION METHODS  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace gamma
    {
        double pdf(double x, double shape, double scale)
        {
            return std::exp(-x / scale) * std::pow(x, shape - 1.) / (std::tgamma(shape) * std::pow(scale, shape));
        }
        double cdf(double x, double shape, double scale)
        {
            return gamma_p(shape, x / scale);
        }
        double logpdf(double x, double shape, double scale)
        {
            return (shape - 1) * std::log(x) - std::lgamma(shape) - shape * std::log(scale) - x / scale;
        }
        double logratio(double xold, double xnew, double shape, double scale)
        {
            return (shape - 1) * std::log(xnew / xold) + (xold - xnew) / scale;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  RAYLEIGH DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace rayleigh
    {
        double pdf(double x, double scale)
        {
            const double scale2 = scale * scale;
            return x * std::exp(-.5 * x * x / scale2) / scale2;
        }
        double cdf(double x, double scale)
        {
            return 1. - std::exp(-.5 * x * x / (scale * scale));
        }
        double logpdf(double x, double scale)
        {
            const double scale2 = scale * scale;
            return std::log(x / scale2) - .5 * x * x / scale2;
        }
        double logratio(double xold, double xnew, double scale)
        {
            return std::log(xnew / xold) + .5 * (xold * xold - xnew * xnew) / (scale * scale);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  WEIBULL DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace weibull
    {
        double _pdf_(double x, double shape, double scale)
        {
            const double z = x / scale;
            return (shape / scale) * std::pow(z, shape - 1.) * std::exp(-std::pow(z, shape));
        }
        double _cdf_(double x, double shape, double scale)
        {
            return 1. - std::exp(-std::pow(x / scale, shape));
        }
        double _logpdf_(double x, double shape, double scale)
        {
            const double z = x / scale;
            return std::log(shape / scale) + (shape - 1.) * std::log(z) - std::pow(z, shape);
        }
        double _logratio_(double xold, double xnew, double shape, double scale)
        {
            return (shape - 1.) * std::log(xnew / xold) + std::pow(xold / scale, shape) - std::pow(xnew / scale, shape);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  EXPONENTIAL DISTRIBUTION METHODS  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace exponential
    {
        double pdf(double x, double rate)
        {
            return rate * std::exp(-rate * x);
        }
        double cdf(double x, double rate)
        {
            return 1. - std::exp(-rate * x);
        }
        double logpdf(double x, double rate)
        {
            return std::log(rate) - rate * x;
        }
        double logratio(double xold, double xnew, double rate)
        {
            return rate * (xold - xnew);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  POISSON DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace poisson
    {
        double pdf(double x, double lambda)
        {
            return std::exp(x * std::log(lambda) - lambda - std::lgamma(x + 1.));
        }
        double cdf(double x, double lambda)
        {
            return 1. - gamma_p(x + 1., lambda);
        }
        double logpdf(double x, double lambda)
        {
            return x * std::log(lambda) - lambda - std::lgamma(x + 1.);
        }
        double logratio(double xold, double xnew, double lambda)
        {
            return (xnew - xold) * std::log(lambda) + std::lgamma(xold + 1.) - std::lgamma(xnew + 1.);
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  BINOMIAL DISTRIBUTION METHODS  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace binomial
    {
        double pdf(double x, double size, double prob)
        {
            return std::exp(log_binom_coeff(size, x) + x * std::log(prob) + (size - x) * std::log(1. - prob));
        }
        double cdf(double x, double size, double prob)
        {
            return ibeta(size - x, x + 1, 1. - prob);
        }
        double logpdf(double x, double size, double prob)
        {
            return log_binom_coeff(size, x) + x * std::log(prob) + (size - x) * std::log(1. - prob);
        }
        double logratio(double xold, double xnew, double size, double prob)
        {
            return logratio_binom_coeff(size, xold, xnew) + (xnew - xold) * std::log(prob / (1. - prob));
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////  NEGATIVE BINOMIAL DISTRIBUTION METHODS  ///////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace negbinomial
    {
        double pdf(double x, double size, double prob)
        {
            return std::exp(log_binom_coeff(x + size - 1, x) + x * std::log(1. - prob) + size * std::log(prob));
        }
        double cdf(double x, double size, double prob)
        {
            return ibeta(size, x + 1, prob);
        }
        double logpdf(double x, double size, double prob)
        {
            return log_binom_coeff(x + size - 1, x) + x * std::log(1. - prob) + size * std::log(prob);
        }
        double logratio(double xold, double xnew, double size, double prob)
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
        double pdf(double x, double proba)
        {
            return proba * std::pow(1. - proba, x);
        }
        double cdf(double x, double proba)
        {
            return 1. - std::pow(1. - proba, x + 1);
        }
        double logpdf(double x, double proba)
        {
            return std::log(proba) + x * std::log(1. - proba);
        }
        double logratio(double xold, double xnew, double proba)
        {
            return (xnew - xold) * std::log(1. - proba);
        }
    }
} // NUMR NAMESPACE
