/*
 * num_recipes.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__NUM_RECIPES__H
#define HOUSEPI__NUM_RECIPES__H

namespace numr
{
    using ushort = unsigned short int;
    using ulong = unsigned long int;

    ulong factorial(ulong n);
    ulong binom_coeff(ushort N, ushort k);
    double log_binom_coeff(ushort N, ushort k);
    double logratio_binom_coeff(ushort N, ushort kold, ushort knew);

    double gamma_p(double shape, double z);
    double lbeta(double a, double b);
    double ibeta(double a, double b, double x);

    namespace normal
    {
        double pdf(double x, double mean, double sd);
        double cdf(double x, double mean, double sd);
        double logpdf(double x, double mean, double sd);
        double logratio(double xold, double xnew, double mean, double sd);
    }

    namespace lognormal
    {
        double pdf(double x, double meanlog, double sdlog);
        double cdf(double x, double meanlog, double sdlog);
        double logpdf(double x, double meanlog, double sdlog);
        double logratio(double xold, double xnew, double meanlog, double sdlog);
    }

    namespace beta
    {
        double pdf(double x, double shape1, double shape2);
        double cdf(double x, double shape1, double shape2);
        double logpdf(double x, double shape1, double shape2);
        double logratio(double xold, double xnew, double shape1, double shape2);
    }

    namespace gamma
    {
        double pdf(double x, double shape, double scale);
        double cdf(double x, double shape, double scale);
        double logpdf(double x, double shape, double scale);
        double logratio(double xold, double xnew, double shape, double scale);
    }

    namespace rayleigh
    {
        double pdf(double x, double scale);
        double cdf(double x, double scale);
        double logpdf(double x, double scale);
        double logratio(double xold, double xnew, double scale);
    }
    
    namespace weibull
    {
        double _pdf_(double x, double shape, double scale);
        double _cdf_(double x, double shape, double scale);
        double _logpdf_(double x, double shape, double scale);
        double _logratio_(double xold, double xnew, double shape, double scale);
    }

    namespace exponential
    {
        double pdf(double x, double rate);
        double cdf(double x, double rate);
        double logpdf(double x, double rate);
        double logratio(double xold, double xnew, double rate);
    }

    namespace poisson
    {
        double pdf(double x, double lambda);
        double cdf(double x, double lambda);
        double logpdf(double x, double lambda);
        double logratio(double xold, double xnew, double lambda);
    }

    namespace binomial
    {
        double pdf(double x, double size, double prob);
        double cdf(double x, double size, double prob);
        double logpdf(double x, double size, double prob);
        double logratio(double xold, double xnew, double size, double prob);
    }

    namespace negbinomial
    {
        double pdf(double x, double size, double prob);
        double cdf(double x, double size, double prob);
        double logpdf(double x, double size, double prob);
        double logratio(double xold, double xnew, double size, double prob);
    }

    namespace geometric
    {
        double pdf(double x, double prob);
        double cdf(double x, double prob);
        double logpdf(double x, double prob);
        double logratio(double xold, double xnew, double prob);
    }
}

#endif // HOUSEPI__NUM_RECIPES__H
