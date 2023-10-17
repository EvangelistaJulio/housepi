// utils.cpp //

#include "utils.h"

using namespace std;

// ////////////////////////////////////////////////////////////////////////// //
// ///////////////////////// AUXILIARY FUNCTIONS //////////////////////////// //
// ////////////////////////////////////////////////////////////////////////// //

/* Cleaner string from forbidden characters*/
void string_cleaner(std::string *token)
{
    for (char r : {'\r', '\n'})
    {
        token->erase(remove(token->begin(), token->end(), r), token->end());
    }
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
                output.emplace(args_stock[i], stod(args_stock[i + 1]));
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
bool stobool(std::string text)
{
    bool output;

    if (text == "TRUE" || text == "true" || text == "True" || text == "1")
    {
        output = true;
    }
    else if (text == "FALSE" || text == "false" || text == "False" || text == "0")
    {
        output = false;
    }
    else
    {
        cerr << "\nERROR: Unrecognized boolean value for: " << text << endl;
        cerr << "Exit...\n"
             << endl;
        exit(EXIT_FAILURE);
    }

    return output;
}

/* function to load delimiter */
char load_delim(std::string text)
{
    char delim;

    if (text.front() == '"' && text.back() == '"')
    {
        // Remove the double quotes and store it as a string
        delim = text.substr(1, text.length() - 2)[0];
    }
    else
    {
        cerr << "\nERROR: Unable to extract delimiter." << endl;
        cerr << "Check: " << text << "\n"
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

///////////////////////////////////////////
// ////////// Uniform_Distr //////////// //
///////////////////////////////////////////

/* uniform probability density function */
double Uniform_Distr::pdf(double x)
{
    (void)x;
    return std::pow(get_arg("up") - get_arg("low"), -1);
}

/* uniform probability density function */
double Uniform_Distr::log_pdf(double x)
{
    (void)x;
    return -std::log(get_arg("up") - get_arg("low"));
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
    return std::exp(-0.5 * std::pow((x - mean) / sd, 2)) * std::pow(sd, -1) * 0.398942280;
}

/* logarithm of the normal probability density function */
double Normal_Distr::log_pdf(double x)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");
    return -0.5 * std::pow((x - mean) / sd, 2) - std::log(sd) - 0.918938533;
}

/* log-ratio of the normal probability density function */
double Normal_Distr::log_ratio(double x, double y)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");
    return 0.5 * (std::pow((x - mean) / sd, 2) - std::pow((y - mean) / sd, 2));
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
    return std::exp(-0.5 * std::pow((std::log(x) - mean) / sd, 2)) * std::pow(x * sd, -1) * 0.398942280;
}

/* logarithm of the lognormal probability density function */
double Lognorm_Distr::log_pdf(double x)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");
    return -0.5 * std::pow((std::log(x) - mean) / sd, 2) - std::log(x) - std::log(sd) - 0.918938533;
}

/* log-ratio of the lognormal probability density function */
double Lognorm_Distr::log_ratio(double x, double y)
{
    double mean = get_arg("mean");
    double sd = get_arg("sd");
    return 0.5 * (std::pow((std::log(x) - mean) / sd, 2) - std::pow((std::log(y) - mean) / sd, 2)) + std::log(x) - std::log(y);
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
    return std::exp(-x / scale) * std::pow(x, shape - 1.) / (std::tgamma(shape) * std::pow(scale, shape));
}

/* logarithm of the gamma probability density function */
double Gamma_Distr::log_pdf(double x)
{
    double shape = get_arg("shape");
    double scale = get_arg("scale");
    return (shape - 1) * std::log(x) - std::lgamma(shape) - shape * std::log(scale) - x / scale;
}

/* log-ratio of the gamma probability density function */
double Gamma_Distr::log_ratio(double x, double y)
{
    double shape = get_arg("shape");
    double scale = get_arg("scale");
    return (shape - 1) * (std::log(y / x)) + (x - y) / scale;
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
    return get_arg("rate") * std::exp(-get_arg("rate") * x);
}

/* logarithm of the exponential probability density function */
double Exponential_Distr::log_pdf(double x)
{
    return std::log(get_arg("rate")) - get_arg("rate") * x;
}

/* log-ratio of the exponential probability density function */
double Exponential_Distr::log_ratio(double x, double y)
{
    return get_arg("rate") * (x - y);
}

/* exponential random generator */
double Exponential_Distr::rand(std::mt19937_64 &rng)
{
    std::exponential_distribution<double> distr_gen(get_arg("rate"));
    return distr_gen(rng);
}
