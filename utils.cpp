// utils.cpp //

#include "utils.h"

using namespace std;

// ////////////////////////////////////////////////////////////////////////// //
// ///////////////////////// AUXILIARY FUNCTIONS //////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

/* Cleaner string from forbidden characters*/
std::string string_cleaner(std::string token)
{
    for (char r : {'\r', '\n'})
    {
        token.erase(remove(token.begin(), token.end(), r), token.end());
    }
    return token;
}

/* Return the dictionary from the string */
kwargs parse_dict_string(std::string args_str, char delim)
{
    std::istringstream args_ss(args_str);
    std::string tmp;
    std::vector<string> args_stock;
    kwargs output = {};

    while (getline(args_ss, tmp, delim))
    {
        args_stock.push_back(tmp);
    }
    if (args_stock.size() == 0)
    {
        output = {}; // empty dictionary
    }
    else if (args_stock.size() % 2 == 1)
    {
        cerr << "\nERROR: Unrecognized dictionary." << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }
    else
    {
        for (size_t i = 0; i < args_stock.size(); i += 2)
        {
            try
            {
                output.push_back({args_stock[i], stod(args_stock[i + 1])});
            }
            catch (const exception &)
            {
                cerr << "\nERROR: Invalid entry of a dictionary." << endl;
                cerr << "Exit...\n"
                     << endl;
                exit(EXIT_FAILURE);
            }
        }
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
        cerr << "\nERROR: Unrecognized boolean value for: " << token << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }
}

/* function to load delimiter */
char load_delim(std::string token)
{
    char delim;

    if (token.front() == '"' && token.back() == '"')
    {
        // Remove the double quotes and store it as a string
        delim = token.substr(1, token.length() - 2)[0];
    }
    else
    {
        cerr << "\nERROR: Unable to extract delimiter." << endl;
        cerr << "Check: " << token << "\n"
             << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    return delim;
}

/* auxiliary function to call for error in build parameters */
void build_column_error(std::string var_name, std::string &line)
{
    cerr << "\nERROR: Unable to extract `" << var_name << "` column." << endl;
    cerr << "Check line: " << line << "\n"
         << endl;
    cerr << "Exit...\n"
         << endl;
    exit(EXIT_FAILURE);
}

/* Linear interpolation for an indexed vector */
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

// ////////////////////////////////////////////////////////////////////////// //
// ////////////////////////// Base Distribution ///////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

// FIXME : this is not foolproof : if someone mixes the order of arguments, the result is not correct

/* Overload constructor for internal arguments */
Base_Distrib::Base_Distrib(kwargs new_args_int)
{
    size_t i = 0;
    args_int = new double[new_args_int.size()]; // dynamic allocation
    for (auto const &x : new_args_int)
    {
        args_int[i] = x.value;
        i++;
    }
    _get_arg_func_ = &Base_Distrib::_get_arg_int_; // assign getter
}

/* Overload constructor for external arguments */
Base_Distrib::Base_Distrib(kwargs_ext new_args_ext)
{
    size_t i = 0;
    args_ext = new const double *[new_args_ext.size()]; // dynamic allocation
    for (auto const &x : new_args_ext)
    {
        // x.first;
        args_ext[i] = x.value;
        i++;
    }
    _get_arg_func_ = &Base_Distrib::_get_arg_ext_; // assign getter
}


/* return informations about the distribution */
std::string Base_Distrib::display(char delim, char dict_delim)
{
    std::string output = get_name();
    for (auto const &x : )
}

///////////////////////////////////////////
// ////////// Uniform_Distr //////////// //
///////////////////////////////////////////

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

////////////////////////////////////////////
// /////////// Normal_Distr ///////////// //
////////////////////////////////////////////

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

////////////////////////////////////////////
// ///////////// LOGNORMAL ////////////// //
////////////////////////////////////////////

/* lognormal probability density function */
double Lognorm_Distr::pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return std::exp(-0.5 * std::pow((std::log(x) - mean) / sd, 2)) * std::pow(x * sd, -1) * 0.398942280;
}

/* logarithm of the lognormal probability density function */
double Lognorm_Distr::log_pdf(double x)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return -0.5 * std::pow((std::log(x) - mean) / sd, 2) - std::log(x) - std::log(sd) - 0.918938533;
}

/* log-ratio of the lognormal probability density function */
double Lognorm_Distr::log_ratio(double x, double y)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    return 0.5 * (std::pow((std::log(x) - mean) / sd, 2) - std::pow((std::log(y) - mean) / sd, 2)) + std::log(x) - std::log(y);
}

/* lognormal random generator */
double Lognorm_Distr::rand(std::mt19937_64 &rng)
{
    double mean = get_arg(0);
    double sd = get_arg(1);
    std::lognormal_distribution<double> distr_gen(mean, sd);
    return distr_gen(rng);
}

////////////////////////////////////////////
// /////////////// GAMMA //////////////// //
////////////////////////////////////////////

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

////////////////////////////////////////////
// //////////// EXPONENTIAL ///////////// //
////////////////////////////////////////////

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
