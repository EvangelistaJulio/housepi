/*
 * distributions.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "distributions.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Incomplete gamma Function  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace auxgammap
{
    /* Lower incomplete gamma function using series expansion. */
    double gamma_cdf_series(double scaled_x, double shape)
    {
        double sum = 1.0 / shape;
        double term = sum;
        for (int i = 1; (term > GAMMAP_EPS * sum) && (i < GAMMAP_ITER); ++i)
        {
            term *= scaled_x / (shape + i);
            sum += term;
        }
        return sum * exp(-scaled_x + shape * log(scaled_x) - lgamma(shape));
    }

    /* Upper incomplete gamma function using continued fraction. */
    double gamma_cdf_cf(double scaled_x, double shape)
    {
        double b = scaled_x + 1.0 - shape;
        double c = 1.0 / GAMMAP_FPMIN;
        double d = 1.0 / b;
        double h = d;
        double delta = 0.0;

        for (int i = 1; ((fabs(delta - 1.0) > GAMMAP_EPS) && (i < GAMMAP_ITER)); ++i)
        {
            double an = -i * (i - shape);
            b += 2.0;
            d = an * d + b;
            if (fabs(d) < GAMMAP_FPMIN)
            {
                d = GAMMAP_FPMIN;
            }
            c = b + an / c;
            if (fabs(c) < GAMMAP_FPMIN)
            {
                c = GAMMAP_FPMIN;
            }
            d = 1.0 / d;
            delta = c * d;
            h *= delta;
        }
        return 1.0 - exp(-scaled_x + shape * log(scaled_x) - lgamma(shape)) * h;
    }

    /* Function to compute the gamma cumulative distribution function. */
    double gamma_cdf(double x, double shape, double scale)
    {
        double scaled_x = x / scale;
        if (scaled_x < shape + 10.)
        {
            return gamma_cdf_series(scaled_x, shape);
        }
        else
        {
            return gamma_cdf_cf(scaled_x, shape);
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Uniform Distribution : methods  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Uniform probability density function. */
double Uniform_Distribution::pdf(double x) const
{
    UNUSED(x);
    return 1.0 / (get_arg(1) - get_arg(0));
}

/* Uniform cumulative density function. */
double Uniform_Distribution::cdf(double x) const
{
    return (x - get_arg(0)) / (get_arg(1) - get_arg(0));
}

/* Uniform probability density function. */
double Uniform_Distribution::logpdf(double x) const
{
    UNUSED(x);
    return -log(get_arg(1) - get_arg(0));
}

/* Uniform probability density function. */
double Uniform_Distribution::logratio(double x, double y) const
{
    UNUSED(x);
    UNUSED(y);
    return 0.0;
}

/* Uniform random generator. */
double Uniform_Distribution::rand(sys_rnd_gen &rng)
{
    std::uniform_real_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Normal probability density function. */
double Normal_Distribution::pdf(double x) const
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return exp(-0.5 * pow((x - mean) / sd, 2)) * 0.398942280 / sd;
}

/* Normal cumulative density function. */
double Normal_Distribution::cdf(double x) const
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.5 * (1. + erf(0.707106781 * (x - mean) / sd));
}

/* Logarithm of the normal probability density function. */
double Normal_Distribution::logpdf(double x) const
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return -0.5 * pow((x - mean) / sd, 2) - log(sd) - 0.918938533;
}

/* log-ratio of the normal probability density function */
double Normal_Distribution::logratio(double x, double y) const
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.5 * (pow((x - mean) / sd, 2) - pow((y - mean) / sd, 2));
}

/* normal random generator */
double Normal_Distribution::rand(sys_rnd_gen &rng)
{
    std::normal_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Lognormal probability density function. */
double Lognorm_Distribution::pdf(double x) const
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return exp(-0.5 * pow((log(x) - meanlog) / sdlog, 2)) * 0.398942280 / (x * sdlog);
}

/* Lognormal probability density function. */
double Lognorm_Distribution::cdf(double x) const
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return 0.5 * (1. + erf(0.707106781 * (log(x) - meanlog) / sdlog));
}

/* logarithm of the lognormal probability density function */
double Lognorm_Distribution::logpdf(double x) const
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return -0.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x) - log(sdlog) - 0.918938533;
}

/* log-ratio of the lognormal probability density function */
double Lognorm_Distribution::logratio(double x, double y) const
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return 0.5 * (pow((log(x) - meanlog) / sdlog, 2) - pow((log(y) - meanlog) / sdlog, 2)) + log(x) - log(y);
}

/* lognormal random generator */
double Lognorm_Distribution::rand(sys_rnd_gen &rng)
{
    std::lognormal_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Gamma : probability density function. */
double Gamma_Distribution::pdf(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
}

/* Gamma : cumulative density function. */
double Gamma_Distribution::cdf(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return auxgammap::gamma_cdf(x, shape, scale);
}

/* Gamma : logarithm of probability density function. */
double Gamma_Distribution::logpdf(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
}

/* Gamma : log-ratio of probability density functions. */
double Gamma_Distribution::logratio(double x, double y) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(y / x)) + (x - y) / scale;
}

/* Gamma : random generator. */
double Gamma_Distribution::rand(sys_rnd_gen &rng)
{
    std::gamma_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Weibull : probability density function. */
double Weibull_Distribution::pdf(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
}

/* Weibull : cumulative density function. */
double Weibull_Distribution::cdf(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return 1. - exp(-pow(x / scale, shape));
}

/* Weibull : logarithm of probability density function. */
double Weibull_Distribution::logpdf(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
}

/* Weibull : log-ratio of probability density function. */
double Weibull_Distribution::logratio(double x, double y) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(y / x)) + pow(x / scale, shape) - pow(y / scale, shape);
}

/* Weibull : random generator */
double Weibull_Distribution::rand(sys_rnd_gen &rng)
{
    std::weibull_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Rayleigh Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Rayleigh : probability density function. */
double Rayleigh_Distribution::pdf(double x) const
{
    double scale = get_arg(0);
    return x * pow(scale, -2) * exp(-0.5 * pow(x / scale, 2));
}

/* Rayleigh : cumulative density function. */
double Rayleigh_Distribution::cdf(double x) const
{
    double scale = get_arg(0);
    return exp(-0.5 * pow(x / scale, 2));
}

/* Rayleigh : logarithm of probability density function. */
double Rayleigh_Distribution::logpdf(double x) const
{
    double scale = get_arg(0);
    return log(x) - 2 * log(scale) - 0.5 * pow(x / scale, 2);
}

/* Rayleigh : log-ratio of probability density function. */
double Rayleigh_Distribution::logratio(double x, double y) const
{
    double scale = get_arg(0);
    return log(y / x) + 0.5 * pow(scale, -2) * (pow(x, 2) - pow(y, 2));
}

/* Rayleigh : random generator. */
double Rayleigh_Distribution::rand(sys_rnd_gen &rng)
{
    std::uniform_real_distribution<double> unif_gen(0, 1);
    return get_arg(0) * sqrt(-2 * unif_gen(rng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Exponential Distribution : methods  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Exponential : probability density function. */
double Exponential_Distribution::pdf(double x) const
{
    double rate = get_arg(0);
    return rate * exp(-rate * x);
}

/* Exponential : cumulative density function. */
double Exponential_Distribution::cdf(double x) const
{
    double rate = get_arg(0);
    return 1. - exp(-rate * x);
}

/* Exponential : logarithm of probability density function. */
double Exponential_Distribution::logpdf(double x) const
{
    double rate = get_arg(0);
    return log(rate) - rate * x;
}

/* Exponential : log-ratio of probability density function. */
double Exponential_Distribution::logratio(double x, double y) const
{
    double rate = get_arg(0);
    return rate * (x - y);
}

/* Exponential : random generator. */
double Exponential_Distribution::rand(sys_rnd_gen &rng)
{
    std::exponential_distribution<double> distr_gen(get_arg(0));
    return distr_gen(rng);
}