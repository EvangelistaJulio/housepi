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
 * numerical_recipes.cpp
 *
 */

#include <string>

#include "numerical_recipes.h"

// namespace containing numerical recipes used to compute distributions
namespace numr
{
    using std::exp;
    using std::fabs;
    using std::lgamma;
    using std::log;
    using std::pow;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  COMBINATORICS FUNCTIONS  ////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    /// Binomial coefficient of N elements in k positions
    unsigned long binom_coeff(unsigned short N, unsigned short k)
    {
        if ((N > 0) && (k != N))
        {
            unsigned long result = N;
            unsigned short stop = (k > N - k) ? k : N - k; // use symmetry
            for (unsigned short i = 1; i < stop; ++i)
                result = (result * (N - i)) / (i + 1);
            return result;
        }

        return 1;
    }

    /// Logarithm of binomial coefficient.
    double log_binom_coeff(unsigned short N, unsigned short k)
    {
        if ((N > 0) && (k != N))
        {
            double result = log(N);
            unsigned short stop = (k > N - k) ? k : N - k; // use symmetry
            for (unsigned short i = 1; i < stop; ++i)
                result += log(N - i) - log(i + 1);
            return result;
        }
        return 0.;
    }

    /// Log-ratio of binomial coefficients.
    double logratio_binom_coeff(unsigned short N, unsigned short ko, unsigned short kn)
    {
        if ((N > 0) && (ko != kn))
        {
            double result = 0.;
            if (ko > kn)
                std::swap(ko, kn);

            for (unsigned short i = ko + 1; i <= kn; ++i)
                result += log(N - i + 1) - log(i);
            return result;
        }
        return 0.;
    }

    /// UPDATE: To better reference this:
    /// Press, Teukolsky, Vetterling, and Flannery
    /// Numerical Recipes 3rd Edition: The Art of Scientific Computing
    namespace gamma_p_impl
    {
        /// Lower regulairzed incomplete gamma function using series expansion.
        double series(double shape, double z)
        {
            double sum = 1. / shape;
            double term = sum;
            for (unsigned m = 1; m <= MAX_ITER; ++m)
            {
                term *= z / (shape + m);
                sum += term;

                if (term < EPSILON * sum)
                    return sum * exp(-z + shape * log(z) - lgamma(shape));
            }

            throw std::overflow_error(
                "Overflow in gamma_p function for x/scale=" + std::to_string(z) +
                ", shape=" + std::to_string(shape) + ".");
        }

        /// Upper regularized incomplete gamma function using continued fraction.
        double continued_fraction(double shape, double z)
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
                if (fabs(d) < MIN_DBL)
                    d = MIN_DBL;
                c = b + an / c;
                if (fabs(c) < MIN_DBL)
                    c = MIN_DBL;
                d = 1. / d;
                delta = c * d;
                h *= delta;

                if (fabs(delta - 1.) < EPSILON)
                    return 1. - exp(-z + shape * log(z) - lgamma(shape)) * h;
            }

            throw std::overflow_error(
                "Overflow in gamma_p function for x/scale=" + std::to_string(z) +
                ", shape=" + std::to_string(shape) + ".");
        }
    }

    /// Regularized incomplete gamma function.
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

    /// UPDATE: to better reference this
    /// (C) Copr. 1986-92 Numerical Recipes Software Vn2-5.
    namespace ibeta_impl
    {
        double continued_fraction(double a, double b, double x)
        {
            const double qab = a + b,
                         qap = a + 1.,
                         qam = a - 1.;
            double c = 1.,
                   d = 1. - qab * x / qap;

            if (fabs(d) < MIN_DBL)
                d = MIN_DBL;
            d = 1. / d;

            double h = d;
            for (unsigned m = 1; m <= MAX_ITER; ++m)
            {
                unsigned m2 = 2 * m;
                double aa = (b - m) * m * x / ((qam + m2) * (a + m2));

                d = 1. + aa * d;
                if (fabs(d) < MIN_DBL)
                    d = MIN_DBL;

                c = 1. + aa / c;
                if (fabs(c) < MIN_DBL)
                    c = MIN_DBL;

                d = 1. / d;
                h *= d * c;
                aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));

                d = 1. + aa * d;
                if (fabs(d) < MIN_DBL)
                    d = MIN_DBL;

                c = 1. + aa / c;
                if (fabs(c) < MIN_DBL)
                    c = MIN_DBL;

                d = 1. / d;

                double del = d * c;
                h *= del;

                if (fabs(del - 1.) < EPSILON)
                    return h;
            }

            throw std::overflow_error(
                "Overflow in ibeta function for x=" + std::to_string(x) +
                ", a=" + std::to_string(a) + ", b=" + std::to_string(b) + ".");
        }
    }

    /// Regularized incomplete beta function.
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
                return exp(lbt + log(ibeta_impl::continued_fraction(a, b, x) / a));
            else
                return 1. - exp(lbt + log(ibeta_impl::continued_fraction(b, a, 1. - x) / b));
        }
    }
} // NUMR NAMESPACE
