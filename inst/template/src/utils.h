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

    void print_stdout(const std::string txt);
    void print_decorator(const std::string title,
                         char decor,
                         const std::string txt = UNDEFINED,
                         bool verbose = true);

    // print progress bar to stdout
    void progress_bar_init(const std::string &txt, bool verbose = true);
    void progress_bar_run(const std::string &txt, int progress, int total, bool verbose = true);
    void progress_bar_end(const std::string &txt, bool verbose = true);

    // print warnings to stderr
    void raise_warning(const std::string &txt, bool &is_raised_once);

    // print error message to stderr
    void raise_error(const std::string &txt);

    std::string percent(int num, int denom);
    std::string truth(bool condition);
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

class Base_Function
{
public:
    //=====Constructors=====//
    Base_Function() {};
    Base_Function(sdargs new_args_int, std::vector<std::string> order = {});
    Base_Function(sdargs_ptr new_args_ptr, std::vector<std::string> order = {});
    virtual ~Base_Function() {};

    // FIXME :
    // need to assign args later in case constructor is empty

    //===== Getter (functional) =====//
    std::string get_name() const { return name; }
    double get_arg(int argID) const { return (this->*_get_arg_func_)(argID); }

    //===== Setter (for internal arguments) =====//
    void set_name(const std::string new_name) { name = new_name; }
    void set_arg(int argID, double new_value) { (this->*_set_arg_func_)(argID, new_value); }

    //===== Show =====//
    std::string show_arg_names() const;
    std::string show_arg_values() const;

private:
    std::string name = UNDEFINED;
    std::vector<std::string> args_name; // ordered list of argument names
    double *args_int = nullptr;         // internal arguments (fixed at initialization)
    const double **args_ptr = nullptr;  // external arguments (pointer to variable)

    double (Base_Function::*_get_arg_func_)(int argID) const;
    double _get_arg_int_(int argID) const { return args_int[argID]; }
    double _get_arg_ptr_(int argID) const { return *(args_ptr[argID]); }

    void (Base_Function::*_set_arg_func_)(int argID, double new_value);
    void _set_arg_int_(int argID, double new_value) { args_int[argID] = new_value; }
    void _set_arg_ptr_(int argID, double new_value)
    {
        UNUSED(new_value);
        screen::raise_error("Attempt to set external argument `" + args_name[argID] + "` within `" + name + "`.");
    }
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
        return INF_DOUBLE;
    }
    virtual double log_pdf(double x)
    {
        UNUSED(x);
        screen::raise_error("Attempt to call a method of an undefined distribution.");
        return INF_DOUBLE;
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
        return INF_DOUBLE;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Distribution : public Base_Distribution
{
public:
    Uniform_Distribution(sdargs new_args = {{"min", (double)0.0}, {"max", (double)1.0}}) : Base_Distribution(new_args, {"min", "max"}) { set_name("uniform"); }
    Uniform_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution(new_args_ptr, {"min", "max"}) { set_name("uniform"); }
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
    Normal_Distribution(sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Base_Distribution(new_args, {"mean", "sd"}) { set_name("normal"); }
    Normal_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution(new_args_ptr, {"mean", "sd"}) { set_name("normal"); }
    ~Normal_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Normal_Distribution_Discrete : public Normal_Distribution
{
public:
    Normal_Distribution_Discrete(sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Normal_Distribution(new_args) { set_name("normal_disc"); }
    Normal_Distribution_Discrete(sdargs_ptr new_args_ptr) : Normal_Distribution(new_args_ptr) { set_name("normal_disc"); }
    ~Normal_Distribution_Discrete() {};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override { return round(Normal_Distribution::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Distribution : public Base_Distribution
{
public:
    Lognorm_Distribution(sdargs new_args = {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}) : Base_Distribution(new_args, {"meanlog", "sdlog"}) { set_name("lognorm"); }
    Lognorm_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution(new_args_ptr, {"meanlog", "sdlog"}) { set_name("lognorm"); }
    ~Lognorm_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Lognorm_Distribution_Discrete : public Lognorm_Distribution
{
public:
    Lognorm_Distribution_Discrete(sdargs new_args = {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}) : Lognorm_Distribution(new_args) { set_name("lognorm_disc"); }
    Lognorm_Distribution_Discrete(sdargs_ptr new_args_ptr) : Lognorm_Distribution(new_args_ptr) { set_name("lognorm_disc"); }
    ~Lognorm_Distribution_Discrete() {};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override { return round(Lognorm_Distribution::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Gamma_Distr : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gamma_Distribution : public Base_Distribution
{
public:
    Gamma_Distribution(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Base_Distribution(new_args, {"shape", "scale"}) { set_name("gamma"); }
    Gamma_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution(new_args_ptr, {"shape", "scale"}) { set_name("gamma"); }
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
    Weibull_Distribution(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Base_Distribution(new_args, {"shape", "scale"}) { set_name("weibull"); }
    Weibull_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution(new_args_ptr, {"shape", "scale"}) { set_name("weibull"); }
    ~Weibull_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Weibull_Distribution_Discrete : public Weibull_Distribution
{
public:
    Weibull_Distribution_Discrete(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Weibull_Distribution(new_args) { set_name("weibull_disc"); }
    Weibull_Distribution_Discrete(sdargs_ptr new_args_ptr) : Weibull_Distribution(new_args_ptr) { set_name("weibull_disc"); }
    ~Weibull_Distribution_Discrete() {};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override { return round(Weibull_Distribution::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Rayleigh_Distr : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Rayleigh_Distribution : public Base_Distribution
{
public:
    Rayleigh_Distribution(sdargs new_args = {{"scale", (double)1.0}}) : Base_Distribution(new_args, {"scale"}) { set_name("rayleigh"); }
    Rayleigh_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution(new_args_ptr, {"scale"}) { set_name("rayleigh"); }
    ~Rayleigh_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Rayleigh_Distribution_Discrete : public Rayleigh_Distribution
{
public:
    Rayleigh_Distribution_Discrete(sdargs new_args = {{"scale", (double)1.0}}) : Rayleigh_Distribution(new_args) { set_name("rayleigh_disc"); }
    Rayleigh_Distribution_Discrete(sdargs_ptr new_args_ptr) : Rayleigh_Distribution(new_args_ptr) { set_name("rayleigh_disc"); }
    ~Rayleigh_Distribution_Discrete() {};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override { return round(Rayleigh_Distribution::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Exponential_Distr : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Exponential_Distribution : public Base_Distribution
{
public:
    Exponential_Distribution(sdargs new_args = {{"rate", (double)1.0}}) : Base_Distribution(new_args, {"rate"}) { set_name("exp"); }
    Exponential_Distribution(sdargs_ptr new_args_ptr) : Base_Distribution(new_args_ptr, {"rate"}) { set_name("exp"); }
    ~Exponential_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Exponential_Distribution_Discrete : public Exponential_Distribution
{
public:
    Exponential_Distribution_Discrete(sdargs new_args = {{"rate", (double)1.0}}) : Exponential_Distribution(new_args) { set_name("exp_disc"); }
    Exponential_Distribution_Discrete(sdargs_ptr new_args_ptr) : Exponential_Distribution(new_args_ptr) { set_name("exp_disc"); }
    ~Exponential_Distribution_Discrete() {};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double logratio(double x, double y) override;
    double rand(sys_rnd_gen &rng) override { return round(Exponential_Distribution::rand(rng)); }
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
template <typename sdargs_int_or_ptr>
Base_Distribution *distribution(const std::string name = UNDEFINED, sdargs_int_or_ptr args = {})
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
    else if (name == "normal_disc")
    {
        output = new Normal_Distribution_Discrete(args);
    }
    else if (name == "lognorm")
    {
        output = new Lognorm_Distribution(args);
    }
    else if (name == "lognorm_disc")
    {
        output = new Lognorm_Distribution_Discrete(args);
    }
    else if (name == "gamma")
    {
        output = new Gamma_Distribution(args);
    }
    else if (name == "weibull")
    {
        output = new Weibull_Distribution(args);
    }
    else if (name == "weibull_disc")
    {
        output = new Weibull_Distribution_Discrete(args);
    }
    else if (name == "rayleigh")
    {
        output = new Rayleigh_Distribution(args);
    }
    else if (name == "rayleigh_disc")
    {
        output = new Rayleigh_Distribution_Discrete(args);
    }
    else if (name == "exp")
    {
        output = new Exponential_Distribution(args);
    }
    else if (name == "exp_disc")
    {
        output = new Exponential_Distribution_Discrete(args);
    }
    else if (name != UNDEFINED)
    {
        screen::raise_error("Unknown distribution `" + name + "`.");
    }

    return output;
}

#endif // UTILS__H