// utils.cpp //

#include "utils.h"

///////////////////////////////////////////////////////////////////////
// /////////////////////// Error Functions ///////////////////////// //
///////////////////////////////////////////////////////////////////////

void exit_with_failure_error(std::string &message)
{
    std::cerr << std::endl;
    std::cerr << "ERROR: " << message << "\nExit...\n"
              << std::endl;
    exit(EXIT_FAILURE);
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////// Base Distribution  ///////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Overload constructor for internal arguments */
Base_Function::Base_Function(sdargs new_args_int, std::vector<std::string> order)
{
    args_name = order;                          // developer assigment of the order of arguments
    args_int = new double[new_args_int.size()]; // dynamic allocation
    for (std::size_t i = 0; i < args_name.size(); i++)
    {
        args_int[i] = new_args_int.at(args_name[i]);
    }
    _get_arg_func_ = &Base_Function::_get_arg_int_; // assign getter
}

/* Overload constructor for external arguments */
Base_Function::Base_Function(sdargs_ext new_args_ext, std::vector<std::string> order)
{
    args_name = order;                                  // developer assigment of the order of arguments
    args_ext = new const double *[new_args_ext.size()]; // dynamic allocation
    for (std::size_t i = 0; i < args_name.size(); i++)
    {
        args_ext[i] = new_args_ext.at(args_name[i]);
    }
    _get_arg_func_ = &Base_Function::_get_arg_ext_; // assign getter
}

/* return informations about the distribution */
std::string Base_Function::display(char delim, char dict_delim)
{
    std::string output = get_name() + delim;
    for (std::size_t i = 0; i < args_name.size() - 1; i++)
    {
        output += args_name[i] + dict_delim + std::to_string(get_arg(i)) + dict_delim;
    }
    std::size_t i = args_name.size() - 1;
    output += args_name[i] + dict_delim + std::to_string(get_arg(i));
    return output;
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////// Uniform Distribution ///////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* uniform probability density function */
double Uniform_Distr::pdf(double x)
{
    double low = get_arg(0);
    double up = get_arg(1);
    (void)x;
    return pow(up - low, -1);
}

/* uniform probability density function */
double Uniform_Distr::log_pdf(double x)
{
    double low = get_arg(0);
    double up = get_arg(1);
    (void)x;
    return -log(up - low);
}

/* uniform probability density function */
double Uniform_Distr::log_ratio(double x, double y)
{
    (void)x;
    (void)y;
    return 0.0;
}

/* uniform random generator */
double Uniform_Distr::rand(std::mt19937_64 &rng)
{
    double low = get_arg(0);
    double up = get_arg(1);
    std::uniform_real_distribution<double> distr_gen(low, up);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////// Normal Distribution ///////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* normal probability density function */
double Normal_Distr::pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return exp(-0.5 * pow((x - mean) / sd, 2)) * pow(sd, -1) * 0.398942280;
}

/* logarithm of the normal probability density function */
double Normal_Distr::log_pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return -0.5 * pow((x - mean) / sd, 2) - log(sd) - 0.918938533;
}

/* log-ratio of the normal probability density function */
double Normal_Distr::log_ratio(double x, double y)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.5 * (pow((x - mean) / sd, 2) - pow((y - mean) / sd, 2));
}

/* normal random generator */
double Normal_Distr::rand(std::mt19937_64 &rng)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    std::normal_distribution<double> distr_gen(mean, sd);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Normal Distribution (discretized) /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

double Normal_Distr_int::pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    int r = round(x);
    return 0.5 * (erf(0.707106781 * (r + 0.5 - mean) / sd) - erf(0.707106781 * (r - 0.5 - mean) / sd));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////// Lognormal Distribution ////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* lognormal probability density function */
double Lognorm_Distr::pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return exp(-0.5 * pow((log(x) - meanlog) / sdlog, 2)) * pow(x * sdlog, -1) * 0.398942280;
}

/* logarithm of the lognormal probability density function */
double Lognorm_Distr::log_pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return -0.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x) - log(sdlog) - 0.918938533;
}

/* log-ratio of the lognormal probability density function */
double Lognorm_Distr::log_ratio(double x, double y)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return 0.5 * (pow((log(x) - meanlog) / sdlog, 2) - pow((log(y) - meanlog) / sdlog, 2)) + log(x) - log(y);
}

/* lognormal random generator */
double Lognorm_Distr::rand(std::mt19937_64 &rng)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    std::lognormal_distribution<double> distr_gen(meanlog, sdlog);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////// Lognormal Distribution (discretized) ////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

double Lognorm_Distr_int::pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    double r = floor(x);
    double aux = (x < 1) ? -1 : erf(0.707106781 * (log(r) - meanlog) / sdlog);
    return 0.5 * (erf(0.707106781 * (log(r + 1) - meanlog) / sdlog) - aux);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Gamma Distribution //////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* gamma probability density function */
double Gamma_Distr::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
}

/* logarithm of the gamma probability density function */
double Gamma_Distr::log_pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
}

/* log-ratio of the gamma probability density function */
double Gamma_Distr::log_ratio(double x, double y)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(y / x)) + (x - y) / scale;
}

/* gamma random generator */
double Gamma_Distr::rand(std::mt19937_64 &rng)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    std::gamma_distribution<double> distr_gen(shape, scale);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Weibull Distribution ////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* weibull probability density function */
double Weibull_Distr::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
}

/* logarithm of the weibull probability density function */
double Weibull_Distr::log_pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
}

/* log-ratio of the weibull probability density function */
double Weibull_Distr::log_ratio(double x, double y)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (log(y / x)) + pow(x / scale, shape) - pow(y / scale, shape);
}

/* weibull random generator */
double Weibull_Distr::rand(std::mt19937_64 &rng)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    std::weibull_distribution<double> distr_gen(shape, scale);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////// Weibull Distribution (discretized) //////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

double Weibull_Distr_int::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    double r = floor(x);
    return exp(-pow(r / scale, shape)) - exp(-pow((r + 1) / scale, shape));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////// Rayleigh Distribution ////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* rayleigh probability density function */
double Rayleigh_Distr::pdf(double x)
{
    double scale = get_arg(0);
    return x * pow(scale, -2) * exp(-0.5 * pow(x / scale, 2));
}

/* logarithm of the rayleigh probability density function */
double Rayleigh_Distr::log_pdf(double x)
{
    double scale = get_arg(0);
    return log(x) - 2 * log(scale) - 0.5 * pow(x / scale, 2);
}

/* log-ratio of the rayleigh probability density function */
double Rayleigh_Distr::log_ratio(double x, double y)
{
    double scale = get_arg(0);
    return log(y / x) + 0.5 * pow(scale, -2) * (pow(x, 2) - pow(y, 2));
}

/* rayleigh random generator */
double Rayleigh_Distr::rand(std::mt19937_64 &rng)
{
    double scale = get_arg(0);
    std::uniform_real_distribution<double> unif_gen(0, 1);
    return scale * sqrt(-2 * unif_gen(rng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////// Rayleigh Distribution (discretized) /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

double Rayleigh_Distr_int::pdf(double x)
{
    double scale = get_arg(0);
    double r = floor(x);
    return exp(-0.5 * pow(r / scale, 2)) - exp(-0.5 * pow((r + 1) / scale, 2));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////// Exponential Distribution //////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* exponential probability density function */
double Exponential_Distr::pdf(double x)
{
    double rate = get_arg(0);
    return rate * exp(-rate * x);
}

/* logarithm of the exponential probability density function */
double Exponential_Distr::log_pdf(double x)
{
    double rate = get_arg(0);
    return log(rate) - rate * x;
}

/* log-ratio of the exponential probability density function */
double Exponential_Distr::log_ratio(double x, double y)
{
    double rate = get_arg(0);
    return rate * (x - y);
}

/* exponential random generator */
double Exponential_Distr::rand(std::mt19937_64 &rng)
{
    double rate = get_arg(0);
    std::exponential_distribution<double> distr_gen(rate);
    return distr_gen(rng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////// Exponential Distribution (discretized) ////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

double Exponential_Distr_int::pdf(double x)
{
    double rate = get_arg(0);
    return exp(-rate * floor(x)) * (1 - exp(-rate));
}

double Exponential_Distr_int::log_pdf(double x)
{
    double rate = get_arg(0);
    return -rate * floor(x) + log(1 - exp(-rate));
}

double Exponential_Distr_int::log_ratio(double x, double y)
{
    double rate = get_arg(0);
    return rate * (floor(x) - floor(y));
}