// user_model.hpp //

#define CONTINUOS_TIME // UPDATE: options to switch to integer time

#ifndef DEF_USER_MODEL__H
#define DEF_USER_MODEL__H

#include <string>
#include <cmath>
#include <random>
#include <unordered_map>
#include <iostream>

// global verbose variable
#ifdef CONTINUOS_TIME
// time values type
typedef double timetype;
#else
typedef int timetype;
// UPDATE : turn distributions into discrete LUT
#endif

const std::string UNDEF = "_none_"; // name for undefined choice

/* Linear interpolation for an indexed array */
double lin_idx_interp(double x_value, const std::vector<double> &y);

///////////////////////////////////////////////////////////////////////
// ////////////////////// Base_Distrib : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

// arguments
typedef std::unordered_map<std::string, double> kwargs;
typedef std::unordered_map<std::string, const double *> kwargs_ext;

class Base_Distrib
{
protected:
    kwargs args = {};         // internal arguments (fixed at initialization)
    kwargs_ext args_ext = {}; // external arguments (pointer to variable)

public:
    //=====Constructors=====//
    Base_Distrib(kwargs new_args)
    {
        for (auto const &x : new_args)
        {
            args.emplace(x.first, x.second);
        }
        _get_arg_func_ = &Base_Distrib::_get_arg_int_;
    }
    Base_Distrib(kwargs_ext new_args_ext)
    {
        for (auto const &x : new_args_ext)
        {
            args_ext.emplace(x.first, x.second);
        }
        _get_arg_func_ = &Base_Distrib::_get_arg_ext_;
    }
    ~Base_Distrib(){};

    // Getter (functional)
    double get_arg(std::string argID) { return (this->*_get_arg_func_)(argID); }

    // Setter (just for fixed)
    // UPDATE: controls on the lists of arguments for given distributions
    void set_arg(std::string argID, double new_value) { args.at(argID) = new_value; }

private:
    double (Base_Distrib::*_get_arg_func_)(std::string argID);
    double _get_arg_int_(std::string argID) { return args.at(argID); }
    double _get_arg_ext_(std::string argID) { return *(args_ext.at(argID)); }
};

///////////////////////////////////////////////////////////////////////
// /////////////////// Distribution_1D : class ///////////////////// //
///////////////////////////////////////////////////////////////////////

class Distribution_1D : public Base_Distrib
{
public:
    using Base_Distrib::Base_Distrib;

    // Virtual functions (dummy)
    virtual double pdf(double x) = 0;
    virtual double log_pdf(double x) = 0;
    virtual double log_ratio(double x, double y) = 0;
    virtual double rand(std::mt19937_64 &rng) = 0;
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
Distribution_1D *distribution_switchboard(std::string name = UNDEF, kwargs_int_or_ext args = {})
{
    if (name == UNDEF)
    {
        return nullptr;
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