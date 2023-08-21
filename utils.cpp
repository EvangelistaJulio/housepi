// utils.cpp //

#include "utils.hpp"

using namespace std;

//----------------------------------------------------------------------
// Linear interpolation for an indexed array
//----------------------------------------------------------------------

double lin_idx_interp(double x_value, const std::vector<double> &y)
{
    size_t x_0 = (int)std::floor(x_value);
    double output = 0.0;

    if (x_0 <= 0)
    { // Saturate to lower value
        output = y[0];
    }
    else if (x_0 >= y.size() - 1)
    { // Saturate to upper value
        output = y.back();
    }
    else
    { // Calculate the interpolated y value
        output = y[x_0] + (x_value - x_0) * (y[x_0 + 1] - y[x_0]);
    }

    return output;
}

/////////////////////////////////////////////////////
//  Distribution : Class
/////////////////////////////////////////////////////

Distribution::Distribution(kwargs new_args)
{
    for (auto const &x : new_args)
    {
        args.emplace(x.first, x.second);
    }
    _get_arg_func_ = &Distribution::_get_arg_int_;
}

Distribution::Distribution(kwargs_ext new_args_ext)
{
    for (auto const &x : new_args_ext)
    {
        args_ext.emplace(x.first, x.second);
    }
    _get_arg_func_ = &Distribution::_get_arg_ext_;
}

///////////////////////////////////////////
// ////////// DISTRIBUTION ///////////// //
///////////////////////////////////////////

/* uniform probability density function */
double Uniform_Distr::pdf(double x)
{
    double low = get_arg("low");
    double up = get_arg("up");
    (void)x;

    // The pdf is defined for low <  up
    // UPDATE : raise a warning for wrong parameters ?
    if (low >= up)
    {
        return 0.0;
    }
    else
    {
        return std::pow(up - low, -1);
    }
}

/* uniform probability density function */
double Uniform_Distr::log_pdf(double x)
{
    double low = get_arg("low");
    double up = get_arg("up");
    (void)x;

    // The pdf is defined for low <  up
    if (low >= up)
    {
        return -INFINITY;
    }
    else
    {
        return -std::log(up - low);
    }
}

/* uniform probability density function */
double Uniform_Distr::log_ratio(double x, double y)
{
    double low = get_arg("low");
    double up = get_arg("up");
    (void)x;
    (void)y;

    // The pdf is defined for low <  up
    if (low >= up)
    {
        return -INFINITY;
    }
    else
    {
        return 0.0;
    }
}

////////////////////////////////////////////
// ////////////// NORMAL //////////////// //
////////////////////////////////////////////

/* normal probability density function */
double Normal_Distr::pdf(double x)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");

    // The pdf is defined for sd > 0.
    if (sd <= 0.0)
    {
        return 0.0;
    }
    else
    {
        return std::exp(-0.5 * std::pow((x - mean) / sd, 2)) * std::pow(sd, -1) * 0.398942280;
    }
}

/* logarithm of the normal probability density function */
double Normal_Distr::log_pdf(double x)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");

    // The logarithm of the pdf is undefined for sd <= 0.
    if (sd <= 0.0)
    {
        return -INFINITY;
    }
    else
    {
        return -0.5 * std::pow((x - mean) / sd, 2) - std::log(sd) - 0.918938533;
    }
}

/* log-ratio of the normal probability density function */
double Normal_Distr::log_ratio(double x, double y)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");

    // The logarithm of the pdf is undefined for sd <= 0.
    if (sd <= 0.0)
    {
        return -INFINITY;
    }
    else
    {
        return 0.5 * (std::pow((y - mean) / sd, 2) - std::pow((x - mean) / sd, 2));
    }
}
////////////////////////////////////////////
// ///////////// LOGNORMAL ////////////// //
////////////////////////////////////////////

/* lognormal probability density function */
double Lognorm_Distr::pdf(double x)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");

    // The pdf is defined for x >= 0 and sd > 0.
    if (x < 0.0 || sd <= 0.0)
    {
        return 0.0;
    }
    else
    {
        return std::exp(-0.5 * std::pow((std::log(x) - mean) / sd, 2)) * std::pow(x * sd, -1) * 0.398942280;
    }
}

/* logarithm of the lognormal probability density function */
double Lognorm_Distr::log_pdf(double x)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");

    // The logarithm of the pdf is undefined for x < 0 and sd <= 0.
    if (x < 0.0 || sd <= 0.0)
    {
        return -INFINITY;
    }
    else
    {
        return -0.5 * std::pow((std::log(x) - mean) / sd, 2) - std::log(x) - std::log(sd) - 0.918938533;
    }
}

/* log-ratio of the lognormal probability density function */
double Lognorm_Distr::log_ratio(double x, double y)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");

    // The logarithm of the pdf is undefined for x < 0 and sd <= 0.
    if (x < 0.0 || sd <= 0.0)
    {
        return -INFINITY;
    }
    else
    {
        return 0.5 * (std::pow((std::log(y) - mean) / sd, 2) - std::pow((std::log(x) - mean) / sd, 2)) + std::log(y) - std::log(x);
    }
}

////////////////////////////////////////////
// /////////////// GAMMA //////////////// //
////////////////////////////////////////////

/* gamma probability density function */
double Gamma_Distr::pdf(double x)
{
    double shape = get_arg("shape");
    double scale = get_arg("scale");

    // The pdf is defined for x >= 0 and shape, scale > 0.
    if (x < 0.0 || shape <= 0.0 || scale <= 0.0)
    {
        return 0.0;
    }
    else
    {
        return std::exp(-x / scale) * std::pow(x, shape - 1.) / (std::tgamma(shape) * std::pow(scale, shape));
    }
}

/* logarithm of the gamma probability density function */
double Gamma_Distr::log_pdf(double x)
{
    double shape = get_arg("shape");
    double scale = get_arg("scale");

    // The logarithm of the pdf is undefined for x < 0 and shape, scale <= 0.
    if (x < 0.0 || shape <= 0.0 || scale <= 0.0)
    {
        return -INFINITY;
    }
    else
    {
        return (shape - 1) * std::log(x) - std::lgamma(shape) - shape * std::log(scale) - x / scale;
    }
}

/* log-ratio of the gamma probability density function */
double Gamma_Distr::log_ratio(double x, double y)
{
    double shape = get_arg("shape");
    double scale = get_arg("scale");

    // The logarithm of the pdf is undefined for x < 0 and shape, scale <= 0.
    if (x < 0.0 || shape <= 0.0 || scale <= 0.0)
    {
        return -INFINITY;
    }
    else
    {
        return (shape - 1) * (std::log(x / y)) + (y - x) / scale;
    }
}

////////////////////////////////////////////
// //////////// EXPONENTIAL ///////////// //
////////////////////////////////////////////

/* exponential probability density function */
double Exponential_Distr::pdf(double x)
{
    double rate = get_arg("rate");

    // The pdf is defined for x >= 0 and rate, scale > 0.
    if (x < 0.0 || rate <= 0.0)
    {
        return 0.0;
    }
    else
    {
        return rate * std::exp(-rate * x);
    }
}

/* logarithm of the exponential probability density function */
double Exponential_Distr::log_pdf(double x)
{
    double rate = get_arg("rate");

    // The pdf is defined for x >= 0 and rate, scale > 0.
    if (x < 0.0 || rate <= 0.0)
    {
        return -INFINITY;
    }
    else
    {
        return std::log(rate) - rate * x;
    }
}

/* log-ratio of the exponential probability density function */
double Exponential_Distr::log_ratio(double x, double y)
{
    double rate = get_arg("rate");

    // The pdf is defined for x >= 0 and rate, scale > 0.
    if (x < 0.0 || rate <= 0.0)
    {
        return -INFINITY;
    }
    else
    {
        return rate * (y - x);
    }
}

//////////////////////////////////////////////////////////
// //////// Time Onset to Time of Transmission //////// //
//////////////////////////////////////////////////////////

/* The Time from Onset of Symptoms to Transmission (TOST) is modeled as a skew-logistic distribution (type I).
This distribution is also a function of the incubation time `t_incub`.
Note: symtoms can occur after infection (pre-symptomatic transmission), in other words `tost` can be negative.
----
Ref:
----
Ferretti, L. et al. The timing of COVID-19 transmission. 2020.09.04.20188516
Preprint at https://doi.org/10.1101/2020.09.04.20188516 (2020).*/

/* Probability distribution function */
double Onset_to_Transm::pdf(timetype tost, timetype t_incub)
{

    double out = 0.0;
    double a = get_arg("alpha");
    double m = get_arg("mu");
    double s = get_arg("sigma");
    double T = get_arg("tau");

    double norm_const = (s / a) * (1. + (t_incub / T - 1.) * std::pow(1 + std::exp(m / s), -a) - t_incub / T * (std::pow(1. + std::exp((m + T) / s), -a)));

    if (tost < 0)
    {
        tost = tost * T / t_incub;
    }
    out = std::exp((m - tost) / s) * std::pow(1. + std::exp((m - tost) / s), -1. - a);

    return out / norm_const;
}

/* Cumulative distribution function */
double Onset_to_Transm::cdf(timetype tost, timetype t_incub)
{
    double out = 0.0;
    double a = get_arg("alpha");
    double m = get_arg("mu");
    double s = get_arg("sigma");
    double T = get_arg("tau");

    double norm_const = (s / a) * (1. + (t_incub / T - 1.) * std::pow(1 + std::exp(m / s), -a) - t_incub / T * (std::pow(1. + std::exp((m + T) / s), -a)));

    if (tost < 0)
    {
        // Integrate from -t_incub to tost
        out = t_incub / T * (std::pow(1. + std::exp((m - tost * T / t_incub) / s), -a) - std::pow(1. + std::exp((m + T) / s), -a));
    }
    else
    {
        // Integrate from -t_incub to 0
        out = t_incub / T * (std::pow(1. + std::exp(m / s), -a) - std::pow(1. + std::exp((m + T) / s), -a));
        // Integrate from 0 to tost
        out += std::pow(1. + std::exp((m - tost) / s), -a) - std::pow(1. + std::exp(m / s), -a);
    }

    return out / norm_const;
}