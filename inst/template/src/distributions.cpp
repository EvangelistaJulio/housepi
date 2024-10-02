/*
 * distributions.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "distributions.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Incomplete gamma Function  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace gammap
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

/* uniform probability density function */
double Uniform_Distribution::pdf(double x)
{
    (void)x;
    return 1.0 / (get_arg(1) - get_arg(0));
}

/* uniform probability density function */
double Uniform_Distribution::log_pdf(double x)
{
    (void)x;
    return -log(get_arg(1) - get_arg(0));
}

/* uniform probability density function */
double Uniform_Distribution::logratio(double x, double y)
{
    (void)x;
    (void)y;
    return 0.0;
}

/* uniform random generator */
double Uniform_Distribution::rand(sys_rnd_gen &rng)
{
    std::uniform_real_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* normal probability density function */
double Normal_Distribution::pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return exp(-0.5 * pow((x - mean) / sd, 2)) * 0.398942280 / sd;
}

/* logarithm of the normal probability density function */
double Normal_Distribution::log_pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return -0.5 * pow((x - mean) / sd, 2) - log(sd) - 0.918938533;
}

/* log-ratio of the normal probability density function */
double Normal_Distribution::logratio(double x, double y)
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
//////////////////////////  Discrete Normal Distribution : methods  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Discretized normal probability density function */
double Discrete_Normal_Distribution::pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    int r = round(x);
    return 0.5 * (erf(0.707106781 * (r + 0.5 - mean) / sd) - erf(0.707106781 * (r - 0.5 - mean) / sd));
}

/* Discretized normal random generator */
double Discrete_Normal_Distribution::rand(sys_rnd_gen &rng)
{
    std::normal_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return round(distr_gen(rng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* lognormal probability density function */
double Lognorm_Distribution::pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return exp(-0.5 * pow((log(x) - meanlog) / sdlog, 2)) * 0.398942280 / (x * sdlog);
}

/* logarithm of the lognormal probability density function */
double Lognorm_Distribution::log_pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return -0.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x) - log(sdlog) - 0.918938533;
}

/* log-ratio of the lognormal probability density function */
double Lognorm_Distribution::logratio(double x, double y)
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
//////////////////////////  Discrete Lognormal Distribution : methods  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Discrete lognormal probability density function. */
double Discrete_Lognorm_Distribution::pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    double r = floor(x);
    double aux = (x < 1) ? -1 : erf(0.707106781 * (log(r) - meanlog) / sdlog);
    return 0.5 * (erf(0.707106781 * (log(r + 1) - meanlog) / sdlog) - aux);
}

/* Discrete lognormal random generator. */
double Discrete_Lognorm_Distribution::rand(sys_rnd_gen &rng)
{
    std::lognormal_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return std::round(distr_gen(rng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* gamma probability density function */
double Gamma_Distribution::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
}

/* logarithm of the gamma probability density function */
double Gamma_Distribution::log_pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
}

/* log-ratio of the gamma probability density function */
double Gamma_Distribution::logratio(double x, double y)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(y / x)) + (x - y) / scale;
}

/* gamma random generator */
double Gamma_Distribution::rand(sys_rnd_gen &rng)
{
    std::gamma_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* weibull probability density function */
double Weibull_Distribution::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
}

/* logarithm of the weibull probability density function */
double Weibull_Distribution::log_pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
}

/* log-ratio of the weibull probability density function */
double Weibull_Distribution::logratio(double x, double y)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(y / x)) + pow(x / scale, shape) - pow(y / scale, shape);
}

/* weibull random generator */
double Weibull_Distribution::rand(sys_rnd_gen &rng)
{
    std::weibull_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Weibull Distribution (discrete) : methods  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Discrete weibull probability density function */
double Discrete_Weibull_Distribution::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    double r = floor(x);
    return exp(-pow(r / scale, shape)) - exp(-pow((r + 1) / scale, shape));
}

/* Discrete weibull random generator */
double Discrete_Weibull_Distribution::rand(sys_rnd_gen &rng)
{
    std::weibull_distribution<double> distr_gen(get_arg(0), get_arg(1));
    return round(distr_gen(rng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Rayleigh Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* rayleigh probability density function */
double Rayleigh_Distribution::pdf(double x)
{
    double scale = get_arg(0);
    return x * pow(scale, -2) * exp(-0.5 * pow(x / scale, 2));
}

/* logarithm of the rayleigh probability density function */
double Rayleigh_Distribution::log_pdf(double x)
{
    double scale = get_arg(0);
    return log(x) - 2 * log(scale) - 0.5 * pow(x / scale, 2);
}

/* log-ratio of the rayleigh probability density function */
double Rayleigh_Distribution::logratio(double x, double y)
{
    double scale = get_arg(0);
    return log(y / x) + 0.5 * pow(scale, -2) * (pow(x, 2) - pow(y, 2));
}

/* rayleigh random generator */
double Rayleigh_Distribution::rand(sys_rnd_gen &rng)
{
    std::uniform_real_distribution<double> unif_gen(0, 1);
    return get_arg(0) * sqrt(-2 * unif_gen(rng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Discrete Rayleigh Distribution : methods  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Discrete rayleigh probability density function */
double Discrete_Rayleigh_Distribution::pdf(double x)
{
    double scale = get_arg(0);
    double r = floor(x);
    return exp(-0.5 * pow(r / scale, 2)) - exp(-0.5 * pow((r + 1) / scale, 2));
}

/* Discrete rayleigh random generator */
double Discrete_Rayleigh_Distribution::rand(sys_rnd_gen &rng)
{
    std::uniform_real_distribution<double> unif_gen(0, 1);
    return round(get_arg(0) * sqrt(-2 * unif_gen(rng)));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Exponential Distribution : methods  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* exponential probability density function */
double Exponential_Distribution::pdf(double x)
{
    double rate = get_arg(0);
    return rate * exp(-rate * x);
}

/* logarithm of the exponential probability density function */
double Exponential_Distribution::log_pdf(double x)
{
    double rate = get_arg(0);
    return log(rate) - rate * x;
}

/* log-ratio of the exponential probability density function */
double Exponential_Distribution::logratio(double x, double y)
{
    double rate = get_arg(0);
    return rate * (x - y);
}

/* exponential random generator */
double Exponential_Distribution::rand(sys_rnd_gen &rng)
{
    std::exponential_distribution<double> distr_gen(get_arg(0));
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Exponential Distribution (discrete) : methods  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Discrete_Exponential_Distribution::pdf(double x)
{
    double rate = get_arg(0);
    return exp(-rate * floor(x)) * (1 - exp(-rate));
}

/* Logarithm of the exponential probability density function. */
double Discrete_Exponential_Distribution::log_pdf(double x)
{
    double rate = get_arg(0);
    return -rate * floor(x) + log(1 - exp(-rate));
}

/* Log-ratio of the discrete exponential probability density function. */
double Discrete_Exponential_Distribution::logratio(double x, double y)
{
    double rate = get_arg(0);
    return rate * (floor(x) - floor(y));
}

/* Discrete exponential random generator. */
double Discrete_Exponential_Distribution::rand(sys_rnd_gen &rng)
{
    std::exponential_distribution<double> distr_gen(get_arg(0));
    return round(distr_gen(rng));
}
