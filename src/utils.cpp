// utils.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "utils.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Stdout Print Functions  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Display a message to stdout (if Verbose). */
void display_stdout(const std::string &txt)
{
    if (Verbose)
    {
        std::cout << txt << std::endl;
    }
}
void display_stdout(const char * txt)
{
    display_stdout(std::string(txt));
}
void display_stdout(const std::ostringstream &oss)
{
    display_stdout(oss.str());
}

/* Display a decorator for a section (if Verbose) */
void display_decorator(const std::string &txt, int skip, char decor)
{
    if (Verbose)
    {
        size_t tot_bar_len = 75;
        size_t half_bar_len = (tot_bar_len- txt.length()) / 2;
        std::ostringstream oss;
        oss << std::string(skip, '\n') << std::string(half_bar_len, decor)
            << txt << std::string(tot_bar_len - txt.length() - half_bar_len, decor);
        std::cout << oss.str() << std::endl;
    }
}
void display_decorator(const char * txt, int skip, char decor)
{
    display_decorator(std::string(txt), skip, decor);
}
void display_decorator(const std::ostringstream &oss, int skip, char decor)
{
    display_decorator(oss.str(), skip, decor);
}

/* Display a percentual progress bar (progress/total). */
void display_progress_bar(const std::string &txt, int progress, int total)
{
    int thin_bar = 100;
    if ((Verbose) && (total >= thin_bar))
    {
        if ((progress % (total / thin_bar)) == 0)
        {
            float percentage = static_cast<float>(progress) / total;
            int bar_width = 50; // - prog_text.length();
            int where = percentage * bar_width;

            std::cerr << txt << ": [" << std::string(where, C_DECOR) << " "
                      << std::string(std::max(0, bar_width - where - 1), ' ') << "] "
                      << std::fixed << std::setprecision(2) << percentage * 100.0 << " %\r";
            std::cerr.flush();
        }
    }
}
void display_progress_bar(const char * txt, int progress, int total)
{
    display_progress_bar(std::string(txt), progress, total);
}
void display_progress_bar(const std::ostringstream &oss, int progress, int total)
{
    display_progress_bar(oss.str(), progress, total);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Warning Functions  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Raise a warning if it wasn't raised before. */
void raise_warning_message(const std::string &txt, bool *is_raised_once)
{
    if (*is_raised_once == false)
    {
        std::cerr << "\nWARNING: " << txt << std::endl;
        *is_raised_once = true;
    }
}
void raise_warning_message(const char * txt, bool *is_raised_once)
{
    raise_warning_message(std::string(txt), is_raised_once);
}
void raise_warning_message(const std::ostringstream &oss, bool *is_raised_once)
{
    raise_warning_message(oss.str(), is_raised_once);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Error Functions  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Display an error message. */
void display_error_message(const std::string &txt)
{
    std::cerr << "\nERROR: " << txt << "\nExit...\n"
              << std::endl;
}
void display_error_message(const char * txt)
{
    display_error_message(std::string(txt));
}
void display_error_message(const std::ostringstream &oss)
{
    display_error_message(oss.str());
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Function : methods  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Overload constructor for INTERNAL arguments. */
Base_Function::Base_Function(sdargs new_args_int, std::vector<std::string> order)
{
    args_name = order;                          // developer assigment of the order of arguments
    args_int = new double[new_args_int.size()]; // dynamic allocation
    for (std::size_t i = 0; i < args_name.size(); i++)
    {
        args_int[i] = new_args_int.at(args_name[i]);
    }
    _get_arg_func_ = &Base_Function::_get_arg_int_; // assign getter
    _set_arg_func_ = &Base_Function::_set_arg_int_; // assign setter
}

/* Overload constructor for EXTERNAL arguments. */
Base_Function::Base_Function(sdargs_ext new_args_ext, std::vector<std::string> order)
{
    args_name = order;                                  // developer assigment of the order of arguments
    args_ext = new const double *[new_args_ext.size()]; // dynamic allocation
    for (std::size_t i = 0; i < args_name.size(); i++)
    {
        args_ext[i] = new_args_ext.at(args_name[i]);
    }
    _get_arg_func_ = &Base_Function::_get_arg_ext_; // assign getter
    _set_arg_func_ = &Base_Function::_set_arg_ext_; // assign setter
}

/* Return information to be displayed about the distribution. */
std::ostringstream Base_Function::message_to_display()
{
    std::ostringstream oss;
    oss << get_name() << C_DELIM;
    for (std::size_t i = 0; i < args_name.size() - 1; i++)
    {
        oss << args_name[i] << C_DICT_DELIM << get_arg(i) << C_DICT_DELIM;
    }
    int i = args_name.size() - 1;
    oss << args_name[i] << C_DICT_DELIM << get_arg(i);
    return oss;
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Uniform Distribution : methods  //////////////////////////////////
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
double Uniform_Distr::logratio(double x, double y)
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
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
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
double Normal_Distr::logratio(double x, double y)
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
/////////////////////////  Normal Distribution (discrete) : methods  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Normal_Distr_int::pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    int r = round(x);
    return 0.5 * (erf(0.707106781 * (r + 0.5 - mean) / sd) - erf(0.707106781 * (r - 0.5 - mean) / sd));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
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
double Lognorm_Distr::logratio(double x, double y)
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
/////////////////////////  Lognormal Distribution (discrete) : methods  /////////////////////////
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
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
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
double Gamma_Distr::logratio(double x, double y)
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
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
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
double Weibull_Distr::logratio(double x, double y)
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
////////////////////////  Weibull Distribution (discrete) : methods  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Weibull_Distr_int::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    double r = floor(x);
    return exp(-pow(r / scale, shape)) - exp(-pow((r + 1) / scale, shape));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Rayleigh Distribution : methods  ///////////////////////////////
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
double Rayleigh_Distr::logratio(double x, double y)
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
////////////////////////  Rayleigh Distribution (discrete) : methods  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Rayleigh_Distr_int::pdf(double x)
{
    double scale = get_arg(0);
    double r = floor(x);
    return exp(-0.5 * pow(r / scale, 2)) - exp(-0.5 * pow((r + 1) / scale, 2));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Exponential Distribution : methods  //////////////////////////////
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
double Exponential_Distr::logratio(double x, double y)
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
////////////////////////  Exponential Distribution (discrete) : methods  ////////////////////////
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

double Exponential_Distr_int::logratio(double x, double y)
{
    double rate = get_arg(0);
    return rate * (floor(x) - floor(y));
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//