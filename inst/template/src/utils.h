/*
 * utils.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef UTILS__H
#define UTILS__H

#include <random>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <chrono>

#ifdef USE_RCPP
#include <Rcpp.h>
#else
#include <iostream>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  32bit/64bit  ///////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Use pseudo-random number generator for 64bit architecture. */
#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
using sys_rnd_gen = std::mt19937_64;
#else
using sys_rnd_gen = std::mt19937;
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////  MACROS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

#define UNUSED(x) (void)(x) // MACRO to avoid unused warnings

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

typedef unsigned short int ushort;
typedef unsigned long int ulong;
typedef std::map<std::string, double> sdargs;             // string to double arguments (internal)
typedef std::map<std::string, const double *> sdargs_ptr; // string to pointer to double arguments (external)
typedef std::map<std::string, std::string> ssargs;        // string to string

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// status aliases
constexpr int SUSCEPTIBLE = 0;     // status of susceptible
constexpr int INFECTIOUS = 1;      // status of infectious
constexpr int MISSING_STATUS = -1; // missing status

// FIXME :NaN and NA to be revised
constexpr double INF_DOUBLE = 1e10;         // numerical infinity value
const std::string UNDEFINED = "None";       // name for undefined choice
const std::string EMPTY_RATE = "__empty__"; // name for Contact rate

// constants
constexpr int MAX_NUMBER_ATTEMPTS = 20;            // max n. of attempts for initialization
constexpr double TARGET_ACCEPT_RATE = 0.234;       // target acceptation rate for MCMC
const std::string IO_MISMATCH = "__io_mismatch__"; // name for undefined choice
constexpr double REFERENCE_VALUE = 1.0;            // auxiliary reference value for the parameters
const std::string REFERENCE = "__reference__";     // name for reference choice
const std::string CONTACT_RATE = "Contact";        // name for Contact rate

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Input/Output Functions  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace screen
{
    constexpr char C_DELIM = ',';      // file delimiter
    constexpr char C_LIST_DELIM = '/'; // dictionary delimiter
    constexpr char C_DECOR = '=';      // decorator for prints

    constexpr int DECOR_WIDTH = 65;   // width of progress bar
    constexpr int PROGBAR_WIDTH = 30; // width of progress bar
    constexpr int PROGBAR_THIN = 17;  // thinning of progress bar

    // deal with undefined or not available value stored as string
    bool is_na(const std::string &txt);

    // output formatter
    std::string percent(int num, int denom);
    std::string truth(bool condition);

    void print_stdout(const std::string &txt);
    void print_decorator(const std::string &title,
                         char decor,
                         const std::string &txt = UNDEFINED,
                         bool verbose = true);

    // print progress bar to stdout
    void progress_bar_init(const std::string &txt, bool verbose = true);
    void progress_bar_run(const std::string &txt, int progress, int total, bool verbose = true);
    void progress_bar_end(const std::string &txt, bool verbose = true);

    // print warnings to stderr
    void raise_warning(const std::string &txt, bool &is_raised_once);

    // print error message to stderr
    void raise_error(const std::string &txt);

}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : class //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME :
// for the moment it is not possible to provide default pointer argumets
// a solution could be point to internal, at that point the whole pointer thing could be collapsed

/* Parent class to all function classes where the arguements are stored in arrays. */
class Base_Function
{
public:
    //=====Constructors=====//
    Base_Function() {};
    Base_Function(const std::string default_name, sdargs default_args_val, std::vector<std::string> default_order = {});
    Base_Function(const std::string default_name, sdargs_ptr default_args_ptr, std::vector<std::string> default_order = {});
    virtual ~Base_Function() {};

    //===== Getter (functional) =====//
    std::string get_name() const { return name; }
    double get_arg(int argID) const { return (this->*_get_arg_method_)(argID); }

    //===== Setter (for internal arguments) =====//
    void set_arg(int argID, double new_value) { (this->*_set_arg_method_)(argID, new_value); }

    //===== Show =====//
    std::string show_arg_names() const;
    std::string show_arg_values() const;

private:
    std::string name = UNDEFINED;
    std::vector<std::string> args_names; // ordered list of argument names
    double *args_val = nullptr;          // internal arguments (fixed at initialization)
    const double **args_ptr = nullptr;   // external arguments (pointer to variable)

    double (Base_Function::*_get_arg_method_)(int) const;
    double _get_arg_val_(int argID) const { return args_val[argID]; }
    double _get_arg_ptr_(int argID) const { return *(args_ptr[argID]); }

    // note : the set_arg is only used for the adaptive
    // for this reason a pointer argument can't be set, but just assigned
    void (Base_Function::*_set_arg_method_)(int, double);
    void _set_arg_val_(int argID, double new_value) { args_val[argID] = new_value; }
    void _set_arg_ptr_(int argID, double new_value)
    {
        UNUSED(new_value);
        std::ostringstream oss;
        oss << "Inconsistent change of pointer argument`" << args_names[argID]
            << "` within `" << get_name() << "`.";
        screen::raise_error(oss.str());
    }

protected:
    //===== Initial non default arguments =====//
    void assign_args(sdargs new_args_val);
    void assign_args(sdargs_ptr new_args_ptr);
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Base_Distribution : public Base_Function
{
public:
    using Base_Function::Base_Function;

    virtual double pdf(double x)
    {
        UNUSED(x);
        screen::raise_error("Attempt to call a method of an undefined distribution.");
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double log_pdf(double x)
    {
        UNUSED(x);
        screen::raise_error("Attempt to call a method of an undefined distribution.");
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double logratio(double x, double y)
    {
        UNUSED(x);
        UNUSED(y);
        return 0.0;
    }
    virtual double rand(sys_rnd_gen &rng)
    {
        UNUSED(rng);
        screen::raise_error("Attempt to call a method of an undefined distribution.");
        return std::numeric_limits<double>::quiet_NaN();
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Distribution : public Base_Distribution
{
public:
    Uniform_Distribution(sdargs new_args = {}) : Base_Distribution("uniform", {{"min", (double)0.0}, {"max", (double)1.0}}, {"min", "max"}) { assign_args(new_args); }
    Uniform_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("uniform", new_args_ptr, {"min", "max"}) {}
    ~Uniform_Distribution() {}

    double pdf(double x) override;
    double log_pdf(double x) override;
    double logratio(double x, double y) override;
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Normal_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Distribution : public Base_Distribution
{
public:
    Normal_Distribution(sdargs new_args = {}) : Base_Distribution("normal", {{"mean", (double)0.0}, {"sd", (double)1.0}}, {"mean", "sd"}) { assign_args(new_args); }
    Normal_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("normal", new_args_ptr, {"mean", "sd"}) {}
    ~Normal_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Normal_Distribution : public Base_Distribution
{
public:
    Discrete_Normal_Distribution(sdargs new_args = {}) : Base_Distribution("normal_discr", {{"mean", (double)0.0}, {"sd", (double)1.0}}, {"mean", "sd"}) { assign_args(new_args); }
    Discrete_Normal_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("normal_discr", new_args_ptr, {"mean", "sd"}) {}
    ~Discrete_Normal_Distribution() {};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Distribution : public Base_Distribution
{
public:
    Lognorm_Distribution(sdargs new_args = {}) : Base_Distribution("lognorm", {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}, {"meanlog", "sdlog"}) { assign_args(new_args); }
    Lognorm_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("lognorm", new_args_ptr, {"meanlog", "sdlog"}) {}
    ~Lognorm_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Lognorm_Distribution : public Base_Distribution
{
public:
    Discrete_Lognorm_Distribution(sdargs new_args = {}) : Base_Distribution("lognorm_discr", {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}, {"meanlog", "sdlog"}) { assign_args(new_args); }
    Discrete_Lognorm_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("lognorm_discr", new_args_ptr, {"meanlog", "sdlog"}) {}
    ~Discrete_Lognorm_Distribution() {};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Gamma_Distr : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gamma_Distribution : public Base_Distribution
{
public:
    Gamma_Distribution(sdargs new_args = {}) : Base_Distribution("gamma", {{"shape", (double)1.0}, {"scale", (double)1.0}}, {"shape", "scale"}) { assign_args(new_args); }
    Gamma_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("gamma", new_args_ptr, {"shape", "scale"}) {}
    ~Gamma_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

// FIXME : the discretization requires boost (for incomplete gamma function)

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Weibull_Distr : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Weibull_Distribution : public Base_Distribution
{
public:
    Weibull_Distribution(sdargs new_args = {}) : Base_Distribution("weibull", {{"shape", (double)1.0}, {"scale", (double)1.0}}, {"shape", "scale"}) { assign_args(new_args); }
    Weibull_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("weibull", new_args_ptr, {"shape", "scale"}) {}
    ~Weibull_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Weibull_Distribution : public Base_Distribution
{
public:
    Discrete_Weibull_Distribution(sdargs new_args = {}) : Base_Distribution("weibull_discr", {{"shape", (double)1.0}, {"scale", (double)1.0}}, {"shape", "scale"}) { assign_args(new_args); }
    Discrete_Weibull_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("weibull_discr", new_args_ptr, {"shape", "scale"}) {}
    ~Discrete_Weibull_Distribution() {};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Rayleigh_Distr : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Rayleigh_Distribution : public Base_Distribution
{
public:
    Rayleigh_Distribution(sdargs new_args = {}) : Base_Distribution("rayleigh", {{"scale", (double)1.0}}, {"scale"}) { assign_args(new_args); }
    Rayleigh_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("rayleigh", new_args_ptr, {"scale"}) {}
    ~Rayleigh_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Rayleigh_Distribution : public Base_Distribution
{
public:
    Discrete_Rayleigh_Distribution(sdargs new_args = {}) : Base_Distribution("rayleigh_discr", {{"scale", (double)1.0}}, {"scale"}) { assign_args(new_args); }
    Discrete_Rayleigh_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("rayleigh_discr", new_args_ptr, {"scale"}) {}
    ~Discrete_Rayleigh_Distribution() {};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Exponential_Distr : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Exponential_Distribution : public Base_Distribution
{
public:
    Exponential_Distribution(sdargs new_args = {}) : Base_Distribution("exp", {{"rate", (double)1.0}}, {"rate"}) { assign_args(new_args); }
    Exponential_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("exp", new_args_ptr, {"rate"}) {}
    ~Exponential_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Exponential_Distribution : public Base_Distribution
{
public:
    Discrete_Exponential_Distribution(sdargs new_args = {}) : Base_Distribution("exp_discr", {{"rate", (double)1.0}}, {"rate"}) { assign_args(new_args); }
    Discrete_Exponential_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution("exp_discr", new_args_ptr, {"rate"}) {}
    ~Discrete_Exponential_Distribution() {};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double logratio(double x, double y) override;
    double rand(sys_rnd_gen &rng) override;
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  distribution switchboard  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return the chosen distribution. */
template <typename sdargs_val_or_ptr>
Base_Distribution *distribution(const std::string name = UNDEFINED, sdargs_val_or_ptr args = {})
{
    Base_Distribution *output = new Base_Distribution();

    if (name == "uniform")
    {
        output = new Uniform_Distribution(args);
    }
    else if (name == "normal")
    {
        output = new Normal_Distribution(args);
    }
    else if (name == "normal_discr")
    {
        output = new Discrete_Normal_Distribution(args);
    }
    else if (name == "lognorm")
    {
        output = new Lognorm_Distribution(args);
    }
    else if (name == "lognorm_discr")
    {
        output = new Discrete_Lognorm_Distribution(args);
    }
    else if (name == "gamma")
    {
        output = new Gamma_Distribution(args);
    }
    else if (name == "weibull")
    {
        output = new Weibull_Distribution(args);
    }
    else if (name == "weibull_discr")
    {
        output = new Discrete_Weibull_Distribution(args);
    }
    else if (name == "rayleigh")
    {
        output = new Rayleigh_Distribution(args);
    }
    else if (name == "rayleigh_discr")
    {
        output = new Discrete_Rayleigh_Distribution(args);
    }
    else if (name == "exp")
    {
        output = new Exponential_Distribution(args);
    }
    else if (name == "exp_discr")
    {
        output = new Discrete_Exponential_Distribution(args);
    }
    else if (!screen::is_na(name))
    {
        screen::raise_error("Unknown distribution `" + name + "`.");
    }

    return output;
}

#endif // UTILS__H