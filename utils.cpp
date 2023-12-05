// utils.cpp //

#include "utils.h"

using namespace std;

///////////////////////////////////////////////////////////////////////
// /////////////////////// Error Functions ///////////////////////// //
///////////////////////////////////////////////////////////////////////

void exit_with_failure_error(std::string &message)
{
    std::cerr << endl;
    std::cerr << "ERROR: " << message << "\nExit...\n"
              << endl;
    exit(EXIT_FAILURE);
}

/* auxiliary function to call for error in build parameters */
void build_column_error(std::string var_name, std::string &line)
{
    std::string err_text = "Unable to extract `" + var_name + "` column.\nCheck line : " + line;
    exit_with_failure_error(err_text);
}

///////////////////////////////////////////////////////////////////////
// ///////////////////// Auxiliary Functions /////////////////////// //
///////////////////////////////////////////////////////////////////////

/* Cleaner string from forbidden characters*/
std::string string_cleaner(std::string token)
{
    for (char r : {'\r', '\n'})
    {
        token.erase(remove(token.begin(), token.end(), r), token.end());
    }
    return token;
}

/* Parse the dictionary from the string */
std::map<std::string, std::string> parse_dict_string(std::string dict_str, char dict_delim)
{
    std::istringstream args_ss(dict_str);
    std::string tmp;
    std::vector<std::string> args_stock;
    std::map<std::string, std::string> output = {};

    // rea the entryes of the dictionary
    while (getline(args_ss, tmp, dict_delim))
    {
        args_stock.push_back(tmp);
    }
    if (args_stock.size() == 0)
    {
        output = {}; // empty dictionary
    }
    else if (args_stock.size() % 2 == 1)
    {
        std::string err_text = "Incomplete dictionary: " + dict_str + ".";
        exit_with_failure_error(err_text);
    }
    else
    {
        for (size_t i = 0; i < args_stock.size(); i += 2)
        {
            try
            {
                output.emplace(args_stock[i], args_stock[i + 1]);
            }
            catch (const exception &)
            {
                std::string err_text = "Invalid entry of dictionary: " + dict_str + ".";
                exit_with_failure_error(err_text);
            }
        }
    }

    return output;
}

/* Parse the scalar arguments dictionary from the string */
kwargs parse_dict_args(std::string args_str, char dict_delim)
{
    std::map<std::string, std::string> input = parse_dict_string(args_str, dict_delim);
    kwargs output = {};
    for (auto const &x : input)
    {
        output.emplace(x.first, stod(x.second));
    }
    return output;
}

/* Parse the option from the string */
std::vector<std::string> parse_option(std::string opt_str, char option_delim)
{
    std::vector<std::string> output = {};
    size_t x = opt_str.find(option_delim);

    // Check if the ":" character was found
    if (x < opt_str.length())
    {
        output.push_back(opt_str.substr(0, x));
        output.push_back(opt_str.substr(x + 1));
    }
    else
    {
        std::string err_text = "Impossible to parse option: " + opt_str + ".";
        exit_with_failure_error(err_text);
    }

    return output;
}

/* function to load boolean variables */
bool stobool(std::string token)
{
    token = string_cleaner(token);

    if (token == "TRUE" || token == "true" || token == "True" || token == "1")
    {
        return true;
    }
    else if (token == "FALSE" || token == "false" || token == "False" || token == "0")
    {
        return false;
    }
    else
    {
        std::string err_text = "Unrecognized boolean value for: " + token + ".";
        exit_with_failure_error(err_text);
        return false;
    }
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
Base_Distrib::Base_Distrib(kwargs new_args_int, std::vector<std::string> order)
{
    args_name = order;                          // developer assigment of the order of arguments
    args_int = new double[new_args_int.size()]; // dynamic allocation
    for (size_t i = 0; i < args_name.size(); i++)
    {
        args_int[i] = new_args_int.at(args_name[i]);
    }
    _get_arg_func_ = &Base_Distrib::_get_arg_int_; // assign getter
}

/* Overload constructor for external arguments */
Base_Distrib::Base_Distrib(kwargs_ext new_args_ext, std::vector<std::string> order)
{
    args_name = order;                                  // developer assigment of the order of arguments
    args_ext = new const double *[new_args_ext.size()]; // dynamic allocation
    for (size_t i = 0; i < args_name.size(); i++)
    {
        args_ext[i] = new_args_ext.at(args_name[i]);
    }
    _get_arg_func_ = &Base_Distrib::_get_arg_ext_; // assign getter
}

/* return informations about the distribution */
std::string Base_Distrib::display(char delim, char dict_delim)
{
    std::string output = get_name() + delim;
    for (size_t i = 0; i < args_name.size() - 1; i++)
    {
        output += args_name[i] + dict_delim + std::to_string(get_arg(i)) + dict_delim;
    }
    size_t i = args_name.size() - 1;
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
    return std::pow(up - low, -1);
}

/* uniform probability density function */
double Uniform_Distr::log_pdf(double x)
{
    double low = get_arg(0);
    double up = get_arg(1);
    (void)x;
    return -std::log(up - low);
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
    return std::exp(-0.5 * std::pow((x - mean) / sd, 2)) * std::pow(sd, -1) * 0.398942280;
}

/* logarithm of the normal probability density function */
double Normal_Distr::log_pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return -0.5 * std::pow((x - mean) / sd, 2) - std::log(sd) - 0.918938533;
}

/* log-ratio of the normal probability density function */
double Normal_Distr::log_ratio(double x, double y)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.5 * (std::pow((x - mean) / sd, 2) - std::pow((y - mean) / sd, 2));
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
    int r = std::round(x);
    return 0.5 * (std::erf(0.707106781 * (r + 0.5 - mean) / sd) - std::erf(0.707106781 * (r - 0.5 - mean) / sd));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////// Lognormal Distribution ////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* lognormal probability density function */
double Lognorm_Distr::pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return std::exp(-0.5 * std::pow((std::log(x) - meanlog) / sdlog, 2)) * std::pow(x * sdlog, -1) * 0.398942280;
}

/* logarithm of the lognormal probability density function */
double Lognorm_Distr::log_pdf(double x)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return -0.5 * std::pow((std::log(x) - meanlog) / sdlog, 2) - std::log(x) - std::log(sdlog) - 0.918938533;
}

/* log-ratio of the lognormal probability density function */
double Lognorm_Distr::log_ratio(double x, double y)
{
    double meanlog = get_arg(0);
    double sdlog = get_arg(1);
    return 0.5 * (std::pow((std::log(x) - meanlog) / sdlog, 2) - std::pow((std::log(y) - meanlog) / sdlog, 2)) + std::log(x) - std::log(y);
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
    double r = std::floor(x);
    double aux = (x < 1) ? -1 : std::erf(0.707106781 * (std::log(r) - meanlog) / sdlog);
    return 0.5 * (std::erf(0.707106781 * (std::log(r + 1) - meanlog) / sdlog) - aux);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Gamma Distribution //////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* gamma probability density function */
double Gamma_Distr::pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return std::exp(-x / scale) * std::pow(x, shape - 1.) / (std::tgamma(shape) * std::pow(scale, shape));
}

/* logarithm of the gamma probability density function */
double Gamma_Distr::log_pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * std::log(x) - std::lgamma(shape) - shape * std::log(scale) - x / scale;
}

/* log-ratio of the gamma probability density function */
double Gamma_Distr::log_ratio(double x, double y)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (std::log(y / x)) + (x - y) / scale;
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
    return (shape / scale) * std::pow((x / scale), shape - 1) * std::exp(-std::pow(x / scale, shape));
}

/* logarithm of the weibull probability density function */
double Weibull_Distr::log_pdf(double x)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return std::log(shape / scale) + (shape - 1) * std::log(x / scale) - std::pow(x / scale, shape);
}

/* log-ratio of the weibull probability density function */
double Weibull_Distr::log_ratio(double x, double y)
{
    double shape = get_arg(0);
    double scale = get_arg(1);
    return (shape - 1) * (std::log(y / x)) + std::pow(x / scale, shape) - std::pow(y / scale, shape);
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
    double r = std::floor(x);
    return std::exp(-std::pow(r / scale, shape)) - std::exp(-std::pow((r + 1) / scale, shape));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////// Rayleigh Distribution ////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* rayleigh probability density function */
double Rayleigh_Distr::pdf(double x)
{
    double scale = get_arg(0);
    return x * std::pow(scale, -2) * std::exp(-0.5 * std::pow(x / scale, 2));
}

/* logarithm of the rayleigh probability density function */
double Rayleigh_Distr::log_pdf(double x)
{
    double scale = get_arg(0);
    return std::log(x) - 2 * std::log(scale) - 0.5 * std::pow(x / scale, 2);
}

/* log-ratio of the rayleigh probability density function */
double Rayleigh_Distr::log_ratio(double x, double y)
{
    double scale = get_arg(0);
    return std::log(y / x) + 0.5 * std::pow(scale, -2) * (std::pow(x, 2) - std::pow(y, 2));
}

/* rayleigh random generator */
double Rayleigh_Distr::rand(std::mt19937_64 &rng)
{
    double scale = get_arg(0);
    std::uniform_real_distribution<double> unif_gen(0, 1);
    return scale * std::sqrt(-2 * unif_gen(rng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////// Rayleigh Distribution (discretized) /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

double Rayleigh_Distr_int::pdf(double x)
{
    double scale = get_arg(0);
    double r = std::floor(x);
    return std::exp(-0.5 * std::pow(r / scale, 2)) - std::exp(-0.5 * std::pow((r + 1) / scale, 2));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////// Exponential Distribution //////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* exponential probability density function */
double Exponential_Distr::pdf(double x)
{
    double rate = get_arg(0);
    return rate * std::exp(-rate * x);
}

/* logarithm of the exponential probability density function */
double Exponential_Distr::log_pdf(double x)
{
    double rate = get_arg(0);
    return std::log(rate) - rate * x;
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
    return std::exp(-rate * std::floor(x)) * (1 - std::exp(-rate));
}

double Exponential_Distr_int::log_pdf(double x)
{
    double rate = get_arg(0);
    return -rate * std::floor(x) + std::log(1 - std::exp(-rate));
}

double Exponential_Distr_int::log_ratio(double x, double y)
{
    double rate = get_arg(0);
    return rate * (std::floor(x) - std::floor(y));
}