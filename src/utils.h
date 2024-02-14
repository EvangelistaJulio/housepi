// utils.h //

#ifndef UTILS__H
#define UTILS__H

//=====Global Macros=====//
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

// NOTE : maybe it's time to give-in and fall back to `boost`

const std::string UNDEFINED = "None";        // name for undefined choice
const std::string REFERENCE = "_reference_"; // name for reference choice

///////////////////////////////////////////////////////////////////////
// //////////////////////////// TYPEDEF //////////////////////////// //
///////////////////////////////////////////////////////////////////////

typedef std::map<std::string, double> sdargs;
typedef std::map<std::string, const double *> sdargs_ext;
typedef std::map<std::string, std::string> ssargs;

//////////////////////////////////////////////////////////////////////////////
// /////////////////////// Raise Error Functions ///////////////////////// //
/////////////////////////////////////////////////////////////////////////////

void exit_with_failure_error(std::string &error_text);

///////////////////////////////////////////////////////////////////////
// ////////////////////// Base_Function : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Base_Function
{
public:
    //=====Constructors=====//
    Base_Function(sdargs new_args_int, std::vector<std::string> order = {});
    Base_Function(sdargs_ext new_args_ext, std::vector<std::string> order = {});
    ~Base_Function(){};

    //=====Getter (functional)=====//
    std::string get_name() const { return name; }
    double get_arg(std::size_t argID) { return (this->*_get_arg_func_)(argID); }

    //=====Setter (for internal arguments)=====//
    void set_name(std::string new_name) { name = new_name; }
    void set_arg(std::size_t argID, double new_value) { args_int[argID] = new_value; }

    std::string display(char delim = DELIM, char dict_delim = DICT_DELIM);

protected:
    std::string name = UNDEFINED;
    std::vector<std::string> args_name;          // ordered list of argument names
    double *args_int = nullptr;        // internal arguments (fixed at initialization)
    const double **args_ext = nullptr; // external arguments (pointer to variable)

private:
    double (Base_Function::*_get_arg_func_)(std::size_t argID);

    // FIXME : gigantic bottleneck // FIXME !! need for translation into array for internal
    double _get_arg_int_(std::size_t argID) { return args_int[argID]; }
    double _get_arg_ext_(std::size_t argID) { return *(args_ext[argID]); }
};

///////////////////////////////////////////////////////////////////////
// /////////////////// Distribution_1D : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Distribution_1D : public Base_Function
{
public:
    using Base_Function::Base_Function;

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
    Undefined_Distr() : Distribution_1D((sdargs){}, {}) {}

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
        std::string err_text = "Attempt to call a method of an undefined distribution.";
        exit_with_failure_error(err_text);
        return INFINITY;
    }
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Uniform_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Uniform_Distr : public Distribution_1D
{
public:
    Uniform_Distr(sdargs new_args = {{"min", (double)0.0}, {"max", (double)1.0}}) : Distribution_1D(new_args, {"min", "max"}) { set_name("uniform"); }
    Uniform_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"min", "max"}) { set_name("uniform"); }
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
    Normal_Distr(sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Distribution_1D(new_args, {"mean", "sd"}) { set_name("normal"); }
    Normal_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"mean", "sd"}) { set_name("normal"); }
    ~Normal_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double log_ratio(double x, double y) override;
    virtual double rand(std::mt19937_64 &rng) override;
};

class Normal_Distr_int : public Normal_Distr
{
public:
    Normal_Distr_int(sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : Normal_Distr(new_args) { set_name("normal_int"); }
    Normal_Distr_int(sdargs_ext new_args_ext) : Normal_Distr(new_args_ext) { set_name("normal_int"); }
    ~Normal_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double log_ratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(std::mt19937_64 &rng) override { return round(Normal_Distr::rand(rng)); }
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Lognorm_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Lognorm_Distr : public Distribution_1D
{
public:
    Lognorm_Distr(sdargs new_args = {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}) : Distribution_1D(new_args, {"meanlog", "sdlog"}) { set_name("lognorm"); }
    Lognorm_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"meanlog", "sdlog"}) { set_name("lognorm"); }
    ~Lognorm_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double log_ratio(double x, double y) override;
    virtual double rand(std::mt19937_64 &rng) override;
};

class Lognorm_Distr_int : public Lognorm_Distr
{
public:
    Lognorm_Distr_int(sdargs new_args = {{"meanlog", (double)0.0}, {"sdlog", (double)1.0}}) : Lognorm_Distr(new_args) { set_name("lognorm_int"); }
    Lognorm_Distr_int(sdargs_ext new_args_ext) : Lognorm_Distr(new_args_ext) { set_name("lognorm_int"); }
    ~Lognorm_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double log_ratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(std::mt19937_64 &rng) override { return round(Lognorm_Distr::rand(rng)); }
};

///////////////////////////////////////////////////////////////////////
// ////////////////////// Gamma_Distr : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

class Gamma_Distr : public Distribution_1D
{
public:
    Gamma_Distr(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Distribution_1D(new_args, {"shape", "scale"}) { set_name("gamma"); }
    Gamma_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"shape", "scale"}) { set_name("gamma"); }
    ~Gamma_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double log_ratio(double x, double y) override;
    virtual double rand(std::mt19937_64 &rng) override;
};

// FIXME : the discretization requies boost (for incomplete gamma function)

///////////////////////////////////////////////////////////////////////
// ///////////////////// Weibull_Distr : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Weibull_Distr : public Distribution_1D
{
public:
    Weibull_Distr(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Distribution_1D(new_args, {"shape", "scale"}) { set_name("weibull"); }
    Weibull_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"shape", "scale"}) { set_name("weibull"); }
    ~Weibull_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double log_ratio(double x, double y) override;
    virtual double rand(std::mt19937_64 &rng) override;
};

class Weibull_Distr_int : public Weibull_Distr
{
public:
    Weibull_Distr_int(sdargs new_args = {{"shape", (double)1.0}, {"scale", (double)1.0}}) : Weibull_Distr(new_args) { set_name("weibull_int"); }
    Weibull_Distr_int(sdargs_ext new_args_ext) : Weibull_Distr(new_args_ext) { set_name("weibull_int"); }
    ~Weibull_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double log_ratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(std::mt19937_64 &rng) override { return round(Weibull_Distr::rand(rng)); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////// Rayleigh_Distr : class ////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

class Rayleigh_Distr : public Distribution_1D
{
public:
    Rayleigh_Distr(sdargs new_args = {{"scale", (double)1.0}}) : Distribution_1D(new_args, {"scale"}) { set_name("rayleigh"); }
    Rayleigh_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"scale"}) { set_name("rayleigh"); }
    ~Rayleigh_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double log_ratio(double x, double y) override;
    virtual double rand(std::mt19937_64 &rng) override;
};

class Rayleigh_Distr_int : public Rayleigh_Distr
{
public:
    Rayleigh_Distr_int(sdargs new_args = {{"scale", (double)1.0}}) : Rayleigh_Distr(new_args) { set_name("rayleigh_int"); }
    Rayleigh_Distr_int(sdargs_ext new_args_ext) : Rayleigh_Distr(new_args_ext) { set_name("rayleigh_int"); }
    ~Rayleigh_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override { return log(pdf(x)); }
    double log_ratio(double x, double y) override { return log(pdf(y) / pdf(x)); }
    double rand(std::mt19937_64 &rng) override { return round(Rayleigh_Distr::rand(rng)); }
};

///////////////////////////////////////////////////////////////////////
// /////////////////// Exponential_Distr : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Exponential_Distr : public Distribution_1D
{
public:
    Exponential_Distr(sdargs new_args = {{"rate", (double)1.0}}) : Distribution_1D(new_args, {"rate"}) { set_name("exp"); }
    Exponential_Distr(sdargs_ext new_args_ext) : Distribution_1D(new_args_ext, {"rate"}) { set_name("exp"); }
    ~Exponential_Distr(){};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double log_ratio(double x, double y) override;
    virtual double rand(std::mt19937_64 &rng) override;
};

class Exponential_Distr_int : public Exponential_Distr
{
public:
    Exponential_Distr_int(sdargs new_args = {{"rate", (double)1.0}}) : Exponential_Distr(new_args) { set_name("exp_int"); }
    Exponential_Distr_int(sdargs_ext new_args_ext) : Exponential_Distr(new_args_ext) { set_name("exp_int"); }
    ~Exponential_Distr_int(){};

    double pdf(double x) override;
    double log_pdf(double x) override;
    double log_ratio(double x, double y) override;
    double rand(std::mt19937_64 &rng) override { return round(Exponential_Distr::rand(rng)); }
};

///////////////////////////////////////////////////////////////////////
// //////////////////// distribution_switchboard /////////////////// //
///////////////////////////////////////////////////////////////////////

/* Return the chosen distribution */
template <typename sdargs_int_or_ext>
Distribution_1D *distribution_switchboard(std::string name = UNDEFINED, sdargs_int_or_ext args = {})
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
    else if (name == "normal_int")
    {
        return new Normal_Distr_int(args);
    }
    else if (name == "lognorm")
    {
        return new Lognorm_Distr(args);
    }
    else if (name == "lognorm_int")
    {
        return new Lognorm_Distr_int(args);
    }
    else if (name == "gamma")
    {
        return new Gamma_Distr(args);
    }
    else if (name == "weibull")
    {
        return new Weibull_Distr(args);
    }
    else if (name == "weibull_int")
    {
        return new Weibull_Distr_int(args);
    }
    else if (name == "rayleigh")
    {
        return new Rayleigh_Distr(args);
    }
    else if (name == "rayleigh_int")
    {
        return new Rayleigh_Distr_int(args);
    }
    else if (name == "exp")
    {
        return new Exponential_Distr(args);
    }
    else if (name == "exp_int")
    {
        return new Exponential_Distr_int(args);
    }
    else
    {
        std::string err_text = "Unknown distribution `" + name + "`.";
        exit_with_failure_error(err_text);
        return new Undefined_Distr();
    }
}

///////////////////////////////////////////////////////////////////////
// //////////////////// Analyzer : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

#endif // UTILS__H