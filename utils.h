// user_model.h //

#ifndef DEF_USER_MODEL__H
#define DEF_USER_MODEL__H


//=====global macros=====//
#define DELIM ','
#define DICT_DELIM '/'
#define OPTION_DELIM ':'
#define UNUSED(x) (void)(x) // MACRO to avoid unused warnings

#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

const std::string UNDEFINED = "None";           // name for undefined choice
const std::string REFERENCE = "_reference_";    // name for reference choice

///////////////////////////////////////////////////////////////////////
// //////////////////////////// TYPEDEF //////////////////////////// //
///////////////////////////////////////////////////////////////////////

typedef std::map<std::string, double> kwargs;
typedef std::map<std::string, const double *> kwargs_ext;

///////////////////////////////////////////////////////////////////////
// /////////////////////// Error Functions ///////////////////////// //
///////////////////////////////////////////////////////////////////////

void exit_with_failure_error(std::string &error_text);
void build_column_error(std::string var_name, std::string &line);

///////////////////////////////////////////////////////////////////////
// ///////////////////// Auxiliary Functions /////////////////////// //
///////////////////////////////////////////////////////////////////////

bool stobool(std::string text);
std::string string_cleaner(std::string token);
std::map<std::string, std::string> parse_dict_string(std::string dict_str, char dict_delim = DICT_DELIM);
std::vector<std::string> parse_option(std::string opt_str, char option_delim = OPTION_DELIM);
kwargs parse_dict_args(std::string args_str, char dict_delim = DICT_DELIM);

///////////////////////////////////////////////////////////////////////
// ////////////////////// Base_Distrib : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Base_Distrib
{
public:
    //=====Constructors=====//
    Base_Distrib(kwargs new_args_int, std::vector<std::string> order = {});
    Base_Distrib(kwargs_ext new_args_ext, std::vector<std::string> order = {});
    ~Base_Distrib(){};

    //=====Getter (functional)=====//
    std::string get_name() const { return name; }
    double get_arg(size_t argID) { return (this->*_get_arg_func_)(argID); }

    //=====Setter (for internal arguments)=====//
    void set_name(std::string new_name) { name = new_name; }
    void set_arg(size_t argID, double new_value) { args_int[argID] = new_value; }

    std::string display(char delim = DELIM, char dict_delim = DICT_DELIM);

protected:
    std::string name = UNDEFINED;
    std::vector<std::string> args_name; // ordered list of argument names
    double *args_int = nullptr;         // internal arguments (fixed at initialization)
    const double **args_ext = nullptr;  // external arguments (pointer to variable)

private:
    double (Base_Distrib::*_get_arg_func_)(size_t argID);

    // FIXME : gigantic bottleneck // FIXME !! need for translation into array for internal
    double _get_arg_int_(size_t argID) { return args_int[argID]; }
    double _get_arg_ext_(size_t argID) { return *(args_ext[argID]); }
};

///////////////////////////////////////////////////////////////////////
// /////////////////// Distribution_1D : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Distribution_1D : public Base_Distrib
{
public:
    using Base_Distrib::Base_Distrib;

    //=====Virtual functions=====//
    virtual double pdf(double x) = 0;
    virtual double log_pdf(double x) = 0;
    virtual double log_ratio(double x, double y) = 0;
    virtual double rand(std::mt19937_64 &rng) = 0;
};

///////////////////////////////////////////////////////////////////////
// //////////////////// Undefined_Distr : class //////////////////// //
///////////////////////////////////////////////////////////////////////

class Undefined_Distr : public Distribution_1D
{
public:
    Undefined_Distr() : Distribution_1D((kwargs){}, {}) {}

    double pdf(double x) override
    {
        UNUSED(x);
        return undefined_call_error();
    }
    double log_pdf(double x) override
    {
        UNUSED(x);
        return undefined_call_error();
    }
    double log_ratio(double x, double y) override
    {
        UNUSED(x);
        UNUSED(y);
        return undefined_call_error();
    }
    double rand(std::mt19937_64 &rng) override
    {
        UNUSED(rng);
        return undefined_call_error();
    }

private:
    double undefined_call_error()
    {
        std::cerr << "\nERROR: Attempt to call a method of an undefined distribution." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
        return INFINITY;
    }
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Uniform_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Uniform_Distr : public Distribution_1D
{
public:
    Uniform_Distr(kwargs new_args = {{"low", (double)0.0}, {"up", (double)1.0}}) : Distribution_1D(new_args, {"low", "up"}) { set_name("uniform"); }
    Uniform_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"low", "up"}) { set_name("uniform"); }
    ~Uniform_Distr() {}

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override;
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Normal_Distr : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

class Normal_Distr : public Distribution_1D
{
public:
    // Constructors
    Normal_Distr(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Distribution_1D(new_args, {"mean", "sd"}) { set_name("normal"); }
    Normal_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"mean", "sd"}) { set_name("normal"); }
    ~Normal_Distr(){};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override;
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Lognorm_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Lognorm_Distr : public Distribution_1D
{
public:
    // Constructors
    Lognorm_Distr(kwargs new_args = {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}) : Distribution_1D(new_args, {"meanlog", "sdlog"}) { set_name("lognorm"); }
    Lognorm_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"meanlog", "sdlog"}) { set_name("lognorm"); }
    ~Lognorm_Distr(){};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override;
};

///////////////////////////////////////////////////////////////////////
// ////////////////////// Gamma_Distr : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

class Gamma_Distr : public Distribution_1D
{
public:
    // Constructors
    Gamma_Distr(kwargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Distribution_1D(new_args, {"shape", "scale"}) { set_name("gamma"); }
    Gamma_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"shape", "scale"}) { set_name("gamma"); }
    ~Gamma_Distr(){};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override;
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Weibull_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Weibull_Distr : public Distribution_1D
{
public:
    // Constructors
    Weibull_Distr(kwargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Distribution_1D(new_args, {"shape", "scale"}) { set_name("weibull"); }
    Weibull_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"shape", "scale"}) { set_name("weibull"); }
    ~Weibull_Distr(){};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override;
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Rayleigh_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Rayleigh_Distr : public Distribution_1D
{
public:
    // Constructors
    Rayleigh_Distr(kwargs new_args = {{"scale", (double)1.0}}) : Distribution_1D(new_args, {"scale"}) { set_name("rayleigh"); }
    Rayleigh_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"scale"}) { set_name("rayleigh"); }
    ~Rayleigh_Distr(){};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override;
};

///////////////////////////////////////////////////////////////////////
// /////////////////// Exponential_Distr : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Exponential_Distr : public Distribution_1D
{
public:
    Exponential_Distr(kwargs new_args = {{"rate", (double)1.0}}) : Distribution_1D(new_args, {"rate"}) { set_name("exp"); }
    Exponential_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"rate"}) { set_name("exp"); }
    ~Exponential_Distr(){};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override;
};

///////////////////////////////////////////////////////////////////////
// //////////////////// distribution_switchboard /////////////////// //
///////////////////////////////////////////////////////////////////////

/* Return the chosen distribution */
template <typename kwargs_int_or_ext>
Distribution_1D *distribution_switchboard(std::string name = UNDEFINED, kwargs_int_or_ext args = {})
{
    if (name == UNDEFINED)
    {
        return new Undefined_Distr();
    }
    else if (name == "uniform")
    {
        return new Uniform_Distr(args);
    }
    else if (name == "normal")
    {
        return new Normal_Distr(args);
    }
    else if (name == "lognorm")
    {
        return new Lognorm_Distr(args);
    }
    else if (name == "gamma")
    {
        return new Gamma_Distr(args);
    }
    else if (name == "weibull")
    {
        return new Weibull_Distr(args);
    }
    else if (name == "rayleigh")
    {
        return new Rayleigh_Distr(args);
    }
    else if (name == "exp")
    {
        return new Exponential_Distr(args);
    }
    else
    {
        std::cerr << "\nERROR: Unknown distribution " << name << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        exit(EXIT_FAILURE);
    }
}

///////////////////////////////////////////////////////////////////////
// //////////////////// Analyzer : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

#endif // DEF_USER_MODEL__H