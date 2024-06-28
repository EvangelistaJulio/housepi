// utils.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

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
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

constexpr int SUSCEPTIBLE = 0; // status of susceptible
constexpr int INFECTIOUS = 1;  // status of infectious

constexpr char C_DELIM = ',';        // file delimiter
constexpr char C_DICT_DELIM = '/';   // dictionary delimiter
constexpr char C_OPTION_DELIM = ':'; // option delimiter
constexpr char C_DECOR = '=';        // decorator for prints

constexpr int DECOR_WIDTH = 65;   // width of progress bar
constexpr int PROGBAR_WIDTH = 30; // width of progress bar
constexpr int PROGBAR_THIN = 241; // thinning of progress bar

constexpr int MAX_NUMBER_ATTEMPTS = 20;      // max n. of attempts for initialization
constexpr double TARGET_ACCEPT_RATE = 0.234; // target acceptation rate for MCMC
constexpr double INF_DOUBLE = 1e10;          // infinity value

constexpr double REFERENCE_VALUE = 1.0;      // auxiliary reference value for the parameters
const std::string UNDEFINED = "None";        // name for undefined choice
const std::string REFERENCE = "_reference_"; // name for reference choice
const std::string CONTACT_RATE = "Contact";  // name for Contact rate

extern bool VERBOSE; // external declaration of verbose global state

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

typedef unsigned short int ushort;
typedef unsigned long int ulong;
typedef std::map<std::string, double> sdargs;             // string to double arguments (internal)
typedef std::map<std::string, const double *> sdargs_ext; // string to pointer to double arguments (external)
typedef std::map<std::string, std::string> ssargs;        // string to string

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Input/Output Functions  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void display_stdout(const std::string &txt);        // base print stdout
void display_stdout(const char *txt);               // overload for char array
void display_stdout(const std::ostringstream &oss); // overload for ostringstream

// print decorator to stdout
void display_decorator(const std::string &txt, int skip = 0, char decor = C_DECOR);
void display_decorator(const char *txt, int skip = 0, char decor = C_DECOR);
void display_decorator(const std::ostringstream &oss, int skip = 0, char decor = C_DECOR);

// print progress bar to stdout
void display_init_progress_bar(const std::string &txt);
void display_progress_bar(const std::string &txt, int progress, int total);
void display_close_progress_bar(const std::string &txt);

// print warnings to stderr
void raise_warning_message(const std::string &txt, bool &is_raised_once);
void raise_warning_message(const char *txt, bool &is_raised_once);
void raise_warning_message(const std::ostringstream &oss, bool &is_raised_once);

// print error message to stderr
void raise_error_message(const std::string &txt);
void raise_error_message(const char *txt);
void raise_error_message(const std::ostringstream &oss);

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : class //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Base_Function
{
public:
    //=====Constructors=====//
    Base_Function(){};
    Base_Function(sdargs new_args_int, std::vector<std::string> order = {});
    Base_Function(sdargs_ext new_args_ext, std::vector<std::string> order = {});
    virtual ~Base_Function(){};

    // FIXME :
    // need to assign args later in case constructor is empty

    //=====Getter (functional)=====//
    std::string get_name() const { return name; }
    double get_arg(int argID) const { return (this->*_get_arg_func_)(argID); }

    //=====Setter (for internal arguments)=====//
    void set_name(std::string new_name) { name = new_name; }
    void set_arg(int argID, double new_value) { (this->*_set_arg_func_)(argID, new_value); }

    std::ostringstream message_to_display();

private:
    std::string name = UNDEFINED;
    std::vector<std::string> args_name; // ordered list of argument names
    double *args_int = nullptr;         // internal arguments (fixed at initialization)
    const double **args_ext = nullptr;  // external arguments (pointer to variable)

    double (Base_Function::*_get_arg_func_)(int argID) const;
    double _get_arg_int_(int argID) const { return args_int[argID]; }
    double _get_arg_ext_(int argID) const { return *(args_ext[argID]); }

    void (Base_Function::*_set_arg_func_)(int argID, double new_value);
    void _set_arg_int_(int argID, double new_value) { args_int[argID] = new_value; }
    void _set_arg_ext_(int argID, double new_value)
    {
        UNUSED(new_value);
        raise_error_message("Attempt to set external argument `" + args_name[argID] + "` within `" + name + "`.");
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Base Distribution_1D : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Distribution_1D : public Base_Function
{
public:
    using Base_Function::Base_Function;

    //=====Virtual functions=====//
    virtual double pdf(double x) = 0;
    virtual double log_pdf(double x) = 0;
    virtual double logratio(double x, double y) = 0;
    virtual double rand(sys_rnd_gen &rng) = 0;
};

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Undefined_Distr : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Undefined_Distr : public Distribution_1D
{
public:
    Undefined_Distr() : Distribution_1D() {}

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
    double logratio(double x, double y) override
    {
        UNUSED(x);
        UNUSED(y);
        return 0.0;
    }
    double rand(sys_rnd_gen &rng) override
    {
        UNUSED(rng);
        return undefined_call_error();
    }

private:
    double undefined_call_error()
    {
        raise_error_message("Attempt to call a method of an undefined distribution.");
        std::exit(EXIT_FAILURE);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Distr : public Distribution_1D
{
public:
    Uniform_Distr(sdargs new_args = {{"min", (double)0.0}, {"max", (double)1.0}}) : Distribution_1D(new_args, {"min", "max"}) { set_name("uniform"); }
    Uniform_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"min", "max"}) { set_name("uniform"); }
    ~Uniform_Distr() {}

    double pdf(double x) override;
    double log_pdf(double x) override;
    double logratio(double x, double y) override;
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Normal_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Distr : public Distribution_1D
{
public:
    Normal_Distr(sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Distribution_1D(new_args, {"mean", "sd"}) { set_name("normal"); }
    Normal_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"mean", "sd"}) { set_name("normal"); }
    ~Normal_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Normal_Distr_int : public Normal_Distr
{
public:
    Normal_Distr_int(sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Normal_Distr(new_args) { set_name("normal_int"); }
    Normal_Distr_int(sdargs_ext new_args_ext) : Normal_Distr(new_args_ext) { set_name("normal_int"); }
    ~Normal_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override { return round(Normal_Distr::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Distr : public Distribution_1D
{
public:
    Lognorm_Distr(sdargs new_args = {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}) : Distribution_1D(new_args, {"meanlog", "sdlog"}) { set_name("lognorm"); }
    Lognorm_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"meanlog", "sdlog"}) { set_name("lognorm"); }
    ~Lognorm_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Lognorm_Distr_int : public Lognorm_Distr
{
public:
    Lognorm_Distr_int(sdargs new_args = {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}) : Lognorm_Distr(new_args) { set_name("lognorm_int"); }
    Lognorm_Distr_int(sdargs_ext new_args_ext) : Lognorm_Distr(new_args_ext) { set_name("lognorm_int"); }
    ~Lognorm_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override { return round(Lognorm_Distr::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Gamma_Distr : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gamma_Distr : public Distribution_1D
{
public:
    Gamma_Distr(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Distribution_1D(new_args, {"shape", "scale"}) { set_name("gamma"); }
    Gamma_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"shape", "scale"}) { set_name("gamma"); }
    ~Gamma_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

// FIXME : the discretization requires boost (for incomplete gamma function)

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Weibull_Distr : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Weibull_Distr : public Distribution_1D
{
public:
    Weibull_Distr(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Distribution_1D(new_args, {"shape", "scale"}) { set_name("weibull"); }
    Weibull_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"shape", "scale"}) { set_name("weibull"); }
    ~Weibull_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Weibull_Distr_int : public Weibull_Distr
{
public:
    Weibull_Distr_int(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Weibull_Distr(new_args) { set_name("weibull_int"); }
    Weibull_Distr_int(sdargs_ext new_args_ext) : Weibull_Distr(new_args_ext) { set_name("weibull_int"); }
    ~Weibull_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override { return round(Weibull_Distr::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Rayleigh_Distr : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Rayleigh_Distr : public Distribution_1D
{
public:
    Rayleigh_Distr(sdargs new_args = {{"scale", (double)1.0}}) : Distribution_1D(new_args, {"scale"}) { set_name("rayleigh"); }
    Rayleigh_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"scale"}) { set_name("rayleigh"); }
    ~Rayleigh_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Rayleigh_Distr_int : public Rayleigh_Distr
{
public:
    Rayleigh_Distr_int(sdargs new_args = {{"scale", (double)1.0}}) : Rayleigh_Distr(new_args) { set_name("rayleigh_int"); }
    Rayleigh_Distr_int(sdargs_ext new_args_ext) : Rayleigh_Distr(new_args_ext) { set_name("rayleigh_int"); }
    ~Rayleigh_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double logratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(sys_rnd_gen &rng) override { return round(Rayleigh_Distr::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Exponential_Distr : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Exponential_Distr : public Distribution_1D
{
public:
    Exponential_Distr(sdargs new_args = {{"rate", (double)1.0}}) : Distribution_1D(new_args, {"rate"}) { set_name("exp"); }
    Exponential_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"rate"}) { set_name("exp"); }
    ~Exponential_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Exponential_Distr_int : public Exponential_Distr
{
public:
    Exponential_Distr_int(sdargs new_args = {{"rate", (double)1.0}}) : Exponential_Distr(new_args) { set_name("exp_int"); }
    Exponential_Distr_int(sdargs_ext new_args_ext) : Exponential_Distr(new_args_ext) { set_name("exp_int"); }
    ~Exponential_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double logratio(double x, double y) override;
    double rand(sys_rnd_gen &rng) override { return round(Exponential_Distr::rand(rng)); }
};

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  distribution_switchboard  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return the chosen distribution. */
template <typename sdargs_int_or_ext>
Distribution_1D *distribution_switchboard(std::string name = UNDEFINED, sdargs_int_or_ext args = {})
{
    Distribution_1D *output = new Undefined_Distr();

    if (name == "uniform")
    {
        output = new Uniform_Distr(args);
    }
    else if (name == "normal")
    {
        output = new Normal_Distr(args);
    }
    else if (name == "normal_int")
    {
        output = new Normal_Distr_int(args);
    }
    else if (name == "lognorm")
    {
        output = new Lognorm_Distr(args);
    }
    else if (name == "lognorm_int")
    {
        output = new Lognorm_Distr_int(args);
    }
    else if (name == "gamma")
    {
        output = new Gamma_Distr(args);
    }
    else if (name == "weibull")
    {
        output = new Weibull_Distr(args);
    }
    else if (name == "weibull_int")
    {
        output = new Weibull_Distr_int(args);
    }
    else if (name == "rayleigh")
    {
        output = new Rayleigh_Distr(args);
    }
    else if (name == "rayleigh_int")
    {
        output = new Rayleigh_Distr_int(args);
    }
    else if (name == "exp")
    {
        output = new Exponential_Distr(args);
    }
    else if (name == "exp_int")
    {
        output = new Exponential_Distr_int(args);
    }
    else if (name != UNDEFINED)
    {
        raise_error_message("Unknown distribution `" + name + "`.");
        std::exit(EXIT_FAILURE);
    }

    return output;
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

#endif // UTILS__H