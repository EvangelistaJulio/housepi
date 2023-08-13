// utils.cpp //

#include "utils.hpp"

using namespace std;

//----------------------------------------------------------------------
// Cumulative probability density functions
//----------------------------------------------------------------------

double tost_normaliz_const(double t_incub, kwargs args)
{
    // Normalization for t_onset_sympt_to_transm_pdf
    double norm_const;
    double a = args.at("alpha");
    double m = args.at("mu");
    double s = args.at("sigma");
    double T = args.at("tau");

    norm_const = 1. + (t_incub / T - 1.) * std::pow(1 + std::exp(m / s), -a) - t_incub / T * (std::pow(1. + std::exp((m + T) / s), -a));
    norm_const *= s / a;
    return norm_const;
}

double t_onset_sympt_to_transm_pdf(double tost, double t_incub, kwargs args)
{
    /* The Time from Onset of Symptoms to Transmission (TOST) is modeled as a skew-logistic distribution (type I).
    This distribution is also a function of the incubation time `t_incub`.
    Note: symtoms can occur after infection (pre-symptomatic transmission), in other words `tost` can be negative.
    ----
    Ref:
    ----
    Ferretti, L. et al. The timing of COVID-19 transmission. 2020.09.04.20188516
    Preprint at https://doi.org/10.1101/2020.09.04.20188516 (2020).*/
    double out = 0.0;
    double a = args.at("alpha");
    double m = args.at("mu");
    double s = args.at("sigma");
    double T = args.at("tau");
    double norm_const = tost_normaliz_const(t_incub, args);

    if (tost < 0)
    {
        tost = tost * T / t_incub;
    }
    out = std::exp((m - tost) / s) * std::pow(1. + std::exp((m - tost) / s), -1. - a);

    return out / norm_const;
}

double t_onset_sympt_to_transm_cdf(double tost, double t_incub, kwargs args)
{
    // Cumulative distribution function for `t_onset_sympt_to_transm_pdf`.
    double out = 0.0;
    double a = args.at("alpha");
    double m = args.at("mu");
    double s = args.at("sigma");
    double T = args.at("tau");
    double norm_const = tost_normaliz_const(t_incub, args);

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

//----------------------------------------------------------------------
// Probability density functions (and their logarithms)
//----------------------------------------------------------------------

double unif_pdf(double x, kwargs args)
{
    // normal probability density function
    double low = args.at("low");
    double up = args.at("up");
    double pdf = 0.0;
    (void)x; // simple fix for Warning

    // The pdf is defined for low <  up
    if (low >= up)
    {
        pdf = 0.0;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        pdf = std::pow(up - low, -1);
    }

    return pdf;
}

double log_unif_pdf(double x, kwargs args)
{
    // normal probability density function
    double low = args.at("low");
    double up = args.at("up");
    double logpdf = 0.0;
    (void)x; // simple fix for Warning

    // The pdf is defined for low <  up
    if (low >= up)
    {
        logpdf = -INFINITY;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        logpdf = -std::log(up - low);
    }

    return logpdf;
}

double normal_pdf(double x, kwargs args)
{
    // normal probability density function
    double mean = args.at("mean");
    double sd = args.at("sd");
    double pdf = 0.0;

    // The pdf is defined for sd > 0.
    if (sd <= 0.0)
    {
        pdf = 0.0;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        pdf = std::exp(-0.5 * std::pow((x - mean) / sd, 2)) * std::pow(sd, -1) * 0.398942280;
    }

    return pdf;
}

double log_normal_pdf(double x, kwargs args)
{
    // logarithm of the normal probability density function
    double mean = args.at("mean");
    double sd = args.at("sd");
    double logpdf = 0.0;

    // The logarithm of the pdf is undefined for sd <= 0.
    if (sd <= 0.0)
    {
        logpdf = -INFINITY;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        logpdf = -0.5 * std::pow((x - mean) / sd, 2) - std::log(sd) - 0.918938533;
    }

    return logpdf;
}

double lognorm_pdf(double x, kwargs args)
{
    // lognormal probability density function
    double mean = args.at("mean");
    double sd = args.at("sd");
    double pdf = 0.0;

    // The pdf is defined for x >= 0 and sd > 0.
    if (x < 0.0 || sd <= 0.0)
    {
        pdf = 0.0;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        pdf = std::exp(-0.5 * std::pow((std::log(x) - mean) / sd, 2)) * std::pow(x * sd, -1) * 0.398942280;
    }

    return pdf;
}

double log_lognorm_pdf(double x, kwargs args)
{
    // logarithm of the lognormal probability density function
    double mean = args.at("mean");
    double sd = args.at("sd");
    double logpdf = 0.0;

    // The logarithm of the pdf is undefined for x < 0 and sd <= 0.
    if (x < 0.0 || sd <= 0.0)
    {
        logpdf = -INFINITY;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        logpdf = -0.5 * std::pow((std::log(x) - mean) / sd, 2) - std::log(x) - std::log(sd) - 0.918938533;
    }

    return logpdf;
}

double gamma_pdf(double x, kwargs args)
{
    // gamma probability density function
    double shape = args.at("shape");
    double scale = args.at("scale");
    double pdf, gamma_value;

    // The pdf is defined for x >= 0 and shape, scale > 0.
    if (x < 0.0 || shape <= 0.0 || scale <= 0.0)
    {
        pdf = 0.0;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        gamma_value = std::tgamma(shape);
        pdf = std::exp(-x / scale) * std::pow(x, shape - 1.) / (gamma_value * std::pow(scale, shape));
    }

    return pdf;
}

double log_gamma_pdf(double x, kwargs args)
{
    // (logarithm of the) gamma probability density function
    double shape = args.at("shape");
    double scale = args.at("scale");
    double logpdf, loggamma_value;

    // The logarithm of the pdf is undefined for x < 0 and shape, scale <= 0.
    if (x < 0.0 || shape <= 0.0 || scale <= 0.0)
    {
        logpdf = -INFINITY;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        loggamma_value = std::lgamma(shape);
        logpdf = loggamma_value - shape * std::log(scale) + (shape - 1) * std::log(x) - x / scale;
    }

    return logpdf;
}

double exp_pdf(double x, kwargs args)
{
    // exponential probability density function
    double rate = args.at("rate");
    double pdf;

    // The pdf is defined for x >= 0 and rate, scale > 0.
    if (x < 0.0 || rate <= 0.0)
    {
        pdf = 0.0;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        pdf = rate * std::exp(-rate * x);
    }

    return pdf;
}

double log_exp_pdf(double x, kwargs args)
{
    // (logarithm of the) exponential probability density function
    double rate = args.at("rate");
    double logpdf;

    // The pdf is defined for x >= 0 and rate, scale > 0.
    if (x < 0.0 || rate <= 0.0)
    {
        logpdf = -INFINITY;
        // NOTE: raise a warning for wrong parameters ?
    }
    else
    {
        logpdf = std::log(rate) - rate * x;
    }

    return logpdf;
}

//----------------------------------------------------------------------
//  Random generators
//----------------------------------------------------------------------

double unif_rand(std::mt19937_64& rand_gen, kwargs args)
{
    // uniform distribution random generator
    double low = args.at("low");
    double up = args.at("up");
    std::uniform_real_distribution<double> distrib(low, up);

    return distrib(rand_gen);
}

double normal_rand(std::mt19937_64& rand_gen, kwargs args)
{
    // normal distribution random generator
    double mean = args.at("mean");
    double sd = args.at("sd");
    std::normal_distribution<double> distrib(mean, sd);

    return distrib(rand_gen);
}

double lognorm_rand(std::mt19937_64& rand_gen, kwargs args)
{
    // lognormal distribution random generator
    double mean = args.at("mean");
    double sd = args.at("sd");
    std::lognormal_distribution<double> distrib(mean, sd);

    return distrib(rand_gen);
}

double gamma_rand(std::mt19937_64& rand_gen, kwargs args)
{
    // gamma distribution random generator
    double shape = args.at("shape");
    double scale = args.at("scale");
    std::gamma_distribution<double> distrib(shape, scale);

    return distrib(rand_gen);
}

double exp_rand(std::mt19937_64& rand_gen, kwargs args)
{
    // exponential distribution random generator
    double rate = args.at("rate");
    std::exponential_distribution<double> distrib(rate);

    return distrib(rand_gen);
}

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