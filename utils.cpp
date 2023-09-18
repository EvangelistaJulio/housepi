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

///////////////////////////////////////////
// ////////// Uniform_Distr //////////// //
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

/* uniform random generator */
double Uniform_Distr::rand(std::mt19937_64 &rng)
{
    std::uniform_real_distribution<double> distr_gen(get_arg("low"), get_arg("up"));
    return distr_gen(rng);
}

////////////////////////////////////////////
// /////////// Normal_Distr ///////////// //
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
        return 0.5 * (std::pow((x - mean) / sd, 2) - std::pow((y - mean) / sd, 2));
    }
}

/* normal random generator */
double Normal_Distr::rand(std::mt19937_64 &rng)
{  
    std::normal_distribution<double> distr_gen(get_arg("mean"), get_arg("sd"));
    return distr_gen(rng);
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
        return 0.5 * (std::pow((std::log(x) - mean) / sd, 2) - std::pow((std::log(y) - mean) / sd, 2)) + std::log(x) - std::log(y);
    }
}

/* lognormal random generator */
double Lognorm_Distr::rand(std::mt19937_64 &rng)
{  
    std::lognormal_distribution<double> distr_gen(get_arg("mean"), get_arg("sd"));
    return distr_gen(rng);
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
        return (shape - 1) * (std::log(y / x)) + (x - y) / scale;
    }
}

/* gamma random generator */
double Gamma_Distr::rand(std::mt19937_64 &rng)
{  
    std::gamma_distribution<double> distr_gen(get_arg("shape"), get_arg("scale"));
    return distr_gen(rng);
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
        return rate * (x - y);
    }
}

/* exponential random generator */
double Exponential_Distr::rand(std::mt19937_64 &rng)
{  
    std::exponential_distribution<double> distr_gen(get_arg("rate"));
    return distr_gen(rng);
}
