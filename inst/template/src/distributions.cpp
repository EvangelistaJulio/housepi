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
/////////////////////////////  Base Distribution : methods  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Internal error for calls to undefined distribution. */
double Base_Distribution::_undefined_call_error_(const std::string &txt, double x, double y) const
{
    UNUSED(x);
    UNUSED(y);
    basics::raise_error("Attempt to call method '" + txt + "' of undefined distribution.");
    return std::numeric_limits<double>::quiet_NaN();
}

/* Error for. */
void Base_Distribution::check_rand_attempts(int &num_attempts) const
{
    if (num_attempts > MAX_NUMBER_ATTEMPTS)
    {
        std::ostringstream oss;
        oss << "Exceeded attempts (N > " << MAX_NUMBER_ATTEMPTS << ")"
            << " in initialization of '" << get_name() << "' "
            << " within range {" << get_range_min() << ", " << get_range_max() << "}.";
        basics::raise_error(oss.str());
    }
    ++num_attempts;
}

/* Truncate distribution at declaration.  */
void Base_Distribution::set_truncation(const dictargs &range)
{
    int trunc_cases = set_range(range);
    if (trunc_cases == 0) // no truncation
    {
        _trunc_method_ = &Base_Distribution::_trunc_dummy_;
        _logtrunc_method_ = &Base_Distribution::_logtrunc_dummy_;
        _shift_method_ = &Base_Distribution::_shift_dummy_;
    }
    else if (trunc_cases == 1) // lower range truncation
    {
        _trunc_method_ = &Base_Distribution::_trunc_min_;
        _logtrunc_method_ = &Base_Distribution::_logtrunc_min_;
        _shift_method_ = &Base_Distribution::_shift_min_;
    }
    else if (trunc_cases == 2) // upper range truncation
    {
        _trunc_method_ = &Base_Distribution::_trunc_max_;
        _logtrunc_method_ = &Base_Distribution::_logtrunc_max_;
        _shift_method_ = &Base_Distribution::_shift_dummy_;
    }
    else // both sides truncation
    {
        _trunc_method_ = &Base_Distribution::_trunc_both_sides_;
        _logtrunc_method_ = &Base_Distribution::_logtrunc_both_sides_;
        _shift_method_ = &Base_Distribution::_shift_dummy_;
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

double Uniform_Distribution::_pdf_(double x) const
{
    UNUSED(x);
    return 1.0;
}

double Uniform_Distribution::_logpdf_(double x) const
{
    UNUSED(x);
    return 0.0;
}

double Uniform_Distribution::_cdf_(double x) const
{
    return x;
}

double Uniform_Distribution::logratio(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.0;
}

double Uniform_Distribution::rand(sys_rnd_gen &prng)
{
    std::uniform_real_distribution<double> distr_gen(get_arg(0), get_arg(1));

    // NOTE : uniform always within range
    return distr_gen(prng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Normal_Distribution::_pdf_(double x) const
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.39894228 * exp(-0.5 * pow((x - mean) / sd, 2)) / sd;
}

double Normal_Distribution::_logpdf_(double x) const
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return -0.5 * pow((x - mean) / sd, 2) - 0.918938533 - log(sd);
}

double Normal_Distribution::_cdf_(double x) const
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.5 + 0.5 * erf(0.707106781 * (x - mean) / sd);
}

double Normal_Distribution::logratio(double xold, double xnew) const
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.5 * (pow((xold - mean) / sd, 2) - pow((xnew - mean) / sd, 2));
}

double Normal_Distribution::rand(sys_rnd_gen &prng)
{
    std::normal_distribution<double> distr_gen(get_arg(0), get_arg(1));
    double new_value = distr_gen(prng);
    int num_attempts = 0;
    while (!valid(new_value))
    {
        check_rand_attempts(num_attempts);
        new_value = distr_gen(prng);
    }
    return new_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Lognorm_Distribution::_pdf_(double x) const
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return exp(-0.5 * pow((log(x) - meanlog) / sdlog, 2)) * 0.39894228 / (x * sdlog);
}

double Lognorm_Distribution::_cdf_(double x) const
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return 0.5 + 0.5 * erf(0.707106781 * (log(x) - meanlog) / sdlog);
}

double Lognorm_Distribution::_logpdf_(double x) const
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return -0.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x * sdlog) - 0.918938533;
}

double Lognorm_Distribution::logratio(double xold, double xnew) const
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return 0.5 * (pow((log(xold) - meanlog) / sdlog, 2) - pow((log(xnew) - meanlog) / sdlog, 2)) + log(xold / xnew);
}

double Lognorm_Distribution::rand(sys_rnd_gen &prng)
{
    std::lognormal_distribution<double> distr_gen(get_arg(0), get_arg(1));
    double new_value = distr_gen(prng);
    int num_attempts = 0;
    while (!valid(new_value))
    {
        check_rand_attempts(num_attempts);
        new_value = distr_gen(prng);
    }
    return new_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Gamma_Distribution::_pdf_(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
}

double Gamma_Distribution::_cdf_(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return auxgammap::gamma_cdf(x, shape, scale);
}

double Gamma_Distribution::_logpdf_(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
}

double Gamma_Distribution::logratio(double xold, double xnew) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(xnew / xold)) + (xold - xnew) / scale;
}

double Gamma_Distribution::rand(sys_rnd_gen &prng)
{
    std::gamma_distribution<double> distr_gen(get_arg(0), get_arg(1));
    double new_value = distr_gen(prng);
    int num_attempts = 0;
    while (!valid(new_value))
    {
        check_rand_attempts(num_attempts);
        new_value = distr_gen(prng);
    }
    return new_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Weibull_Distribution::_pdf_(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
}

double Weibull_Distribution::_cdf_(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return 1. - exp(-pow(x / scale, shape));
}

double Weibull_Distribution::_logpdf_(double x) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
}

double Weibull_Distribution::logratio(double xold, double xnew) const
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1.) * (log(xnew / xold)) + pow(xold / scale, shape) - pow(xnew / scale, shape);
}

double Weibull_Distribution::rand(sys_rnd_gen &prng)
{
    std::weibull_distribution<double> distr_gen(get_arg(0), get_arg(1));
    double new_value = distr_gen(prng);
    int num_attempts = 0;
    while (!valid(new_value))
    {
        check_rand_attempts(num_attempts);
        new_value = distr_gen(prng);
    }
    return new_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Rayleigh Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Rayleigh : probability density function. */
double Rayleigh_Distribution::_pdf_(double x) const
{
    double scale = get_arg(0);
    return x * pow(scale, -2) * exp(-0.5 * pow(x / scale, 2));
}

/* Rayleigh : cumulative density function. */
double Rayleigh_Distribution::_cdf_(double x) const
{
    double scale = get_arg(0);
    return exp(-0.5 * pow(x / scale, 2));
}

/* Rayleigh : logarithm of probability density function. */
double Rayleigh_Distribution::_logpdf_(double x) const
{
    double scale = get_arg(0);
    return log(x) - 2 * log(scale) - 0.5 * pow(x / scale, 2);
}

/* Rayleigh : log-ratio of probability density function. */
double Rayleigh_Distribution::logratio(double xold, double xnew) const
{
    double scale = get_arg(0);
    return log(xnew / xold) + 0.5 * pow(scale, -2) * (pow(xold, 2) - pow(xnew, 2));
}

/* Rayleigh : random generator. */
double Rayleigh_Distribution::rand(sys_rnd_gen &prng)
{
    std::uniform_real_distribution<double> unif_gen(0, 1);
    double new_value = get_arg(0) * sqrt(-2 * unif_gen(prng));
    int num_attempts = 0;
    while (!valid(new_value))
    {
        check_rand_attempts(num_attempts);
        new_value = get_arg(0) * sqrt(-2 * unif_gen(prng));
    }
    return new_value;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Exponential Distribution : methods  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Exponential : probability density function. */
double Exponential_Distribution::_pdf_(double x) const
{
    double rate = get_arg(0);
    return rate * exp(-rate * x);
}

/* Exponential : cumulative density function. */
double Exponential_Distribution::_cdf_(double x) const
{
    double rate = get_arg(0);
    return 1. - exp(-rate * x);
}

/* Exponential : logarithm of probability density function. */
double Exponential_Distribution::_logpdf_(double x) const
{
    double rate = get_arg(0);
    return log(rate) - rate * x;
}

/* Exponential : log-ratio of probability density function. */
double Exponential_Distribution::logratio(double xold, double xnew) const
{
    double rate = get_arg(0);
    return rate * (xold - xnew);
}

/* Exponential : random generator. */
double Exponential_Distribution::rand(sys_rnd_gen &prng)
{
    std::exponential_distribution<double> distr_gen(get_arg(0));
    double new_value = distr_gen(prng);
    int num_attempts = 0;
    while (!valid(new_value))
    {
        check_rand_attempts(num_attempts);
        new_value = distr_gen(prng);
    }
    return new_value;
}