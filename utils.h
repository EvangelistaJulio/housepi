// user_model.h //

#ifndef DEF_USER_MODEL__H
#define DEF_USER_MODEL__H

#define CONTINUOS_TIME      // UPDATE: options to switch to integer time
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

const std::string UNDEFINED = "_none_"; // name for undefined choice

///////////////////////////////////////////////////////////////////////
// //////////////////////////// TYPEDEF //////////////////////////// //
///////////////////////////////////////////////////////////////////////

#ifdef CONTINUOS_TIME
typedef double timetype; // time values type
#else
typedef int timetype;
// UPDATE : turn distributions into discrete LUT
#endif

struct key_value {
    std::string key;
    double value;
};

struct key_address {
    std::string key;
    const double * value;
};

typedef std::vector<key_value> kwargs;
typedef std::vector<key_address> kwargs_ext;

///////////////////////////////////////////////////////////////////////
// ///////////////////// Auxiliary Functions /////////////////////// //
///////////////////////////////////////////////////////////////////////

std::string string_cleaner(std::string token);
kwargs parse_dict_string(std::string args_str, char delim);
bool stobool(std::string text);
char load_delim(std::string text);
void build_column_error(std::string var_name, std::string &line);
double lin_idx_interp(double x_value, const std::vector<double> &y);

///////////////////////////////////////////////////////////////////////
// ////////////////////// Base_Distrib : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Base_Distrib
{
protected:
    double * args_int;   // internal arguments (fixed at initialization)
    const double ** args_ext ; // external arguments (pointer to variable)

public:
    //=====Constructors=====//
    Base_Distrib(kwargs new_args_int);      // FIXME : current definition does not work because of lost order
    Base_Distrib(kwargs_ext new_args_ext);
    ~Base_Distrib(){};

    //=====Getter (functional)=====//
    double get_arg(size_t argID) { return (this->*_get_arg_func_)(argID); }

    //=====Setter (for internal arguments)=====//
    void set_arg(size_t argID, double new_value) { args_int[argID] = new_value; }

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
private:
    double undefined_call_error()
    {
        std::cerr << "\nERROR: Attempt to call a method of an undefined distribution." << std::endl;
        std::cerr << "Exit...\n"
                  << std::endl;
        std::exit(EXIT_FAILURE);
        return INFINITY;
    }

public:
    Undefined_Distr() : Distribution_1D((kwargs){}) {}

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
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Uniform_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Uniform_Distr : public Distribution_1D
{
public:
    Uniform_Distr(kwargs new_args = {{"low", (double)0.0}, {"up", (double)1.0}}) : Distribution_1D(new_args) {}
    Uniform_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext) {}
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
    Normal_Distr(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Distribution_1D(new_args) {}
    Normal_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext) {}
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
    Lognorm_Distr(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Distribution_1D(new_args) {}
    Lognorm_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext) {}
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
    Gamma_Distr(kwargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Distribution_1D(new_args) {}
    Gamma_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext) {}
    ~Gamma_Distr(){};

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
    Exponential_Distr(kwargs new_args = {{"rate", (double)1.0}}) : Distribution_1D(new_args) {}
    Exponential_Distr(kwargs_ext new_args_ext) : Distribution_1D(new_args_ext) {}
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

#endif // DEF_USER_MODEL__H