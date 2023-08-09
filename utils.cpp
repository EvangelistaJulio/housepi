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
    double a = args["alpha"];
    double m = args["mu"];
    double s = args["sigma"];
    double T = args["tau"];

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
    double a = args["alpha"];
    double m = args["mu"];
    double s = args["sigma"];
    double T = args["tau"];
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
    double a = args["alpha"];
    double m = args["mu"];
    double s = args["sigma"];
    double T = args["tau"];
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
    double low = args["low"];
    double up = args["up"];
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
    double low = args["low"];
    double up = args["up"];
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
    double mean = args["mean"];
    double sd = args["sd"];
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
    double mean = args["mean"];
    double sd = args["sd"];
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
    double mean = args["mean"];
    double sd = args["sd"];
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
    double mean = args["mean"];
    double sd = args["sd"];
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
    double shape = args["shape"];
    double scale = args["scale"];
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

//----------------------------------------------------------------------
//  Random generators
//----------------------------------------------------------------------

double unif_rand(std::mt19937_64 &rand_gen, kwargs args)
{
    // uniform distribution random generator
    double low = args["low"];
    double up = args["up"];
    std::uniform_real_distribution<double> distrib(low, up);

    return distrib(rand_gen);
}

double normal_rand(std::mt19937_64 &rand_gen, kwargs args)
{
    // normal distribution random generator
    double mean = args["mean"];
    double sd = args["sd"];
    std::normal_distribution<double> distrib(mean, sd);

    return distrib(rand_gen);
}

double lognorm_rand(std::mt19937_64 &rand_gen, kwargs args)
{
    // lognormal distribution random generator
    double mean = args["mean"];
    double sd = args["sd"];
    std::lognormal_distribution<double> distrib(mean, sd);

    return distrib(rand_gen);
}

double gamma_rand(std::mt19937_64 &rand_gen, kwargs args)
{
    // gamma distribution random generator
    double shape = args["shape"];
    double scale = args["scale"];
    std::gamma_distribution<double> distrib(shape, scale);

    return distrib(rand_gen);
}

double exp_rand(std::mt19937_64 &rand_gen, kwargs args)
{
    // exponential distribution random generator
    double rate = args["rate"];
    std::exponential_distribution<double> distrib(rate);

    return distrib(rand_gen);
}