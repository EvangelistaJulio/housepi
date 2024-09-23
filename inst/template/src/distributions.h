/*
 * distributions.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef DISTRIBUTIONS__H
#define DISTRIBUTIONS__H

#include <random>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "utils.h"

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
    Uniform_Distribution(dictargs new_args = {}) : Base_Distribution("uniform", {{"min", 0.0}, {"max", 1.0}}, {"min", "max"}) { assign_args(new_args); }
    Uniform_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("uniform", new_args_ptr, {"min", "max"}) {}
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
    Normal_Distribution(dictargs new_args = {}) : Base_Distribution("normal", {{"mean", 0.0}, {"sd", 1.0}}, {"mean", "sd"}) { assign_args(new_args); }
    Normal_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("normal", new_args_ptr, {"mean", "sd"}) {}
    ~Normal_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Normal_Distribution : public Base_Distribution
{
public:
    Discrete_Normal_Distribution(dictargs new_args = {}) : Base_Distribution("normal_discr", {{"mean", 0.0}, {"sd", 1.0}}, {"mean", "sd"}) { assign_args(new_args); }
    Discrete_Normal_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("normal_discr", new_args_ptr, {"mean", "sd"}) { set_range({{"min", 0}}); }
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
    Lognorm_Distribution(dictargs new_args = {}) : Base_Distribution("lognorm", {{"meanlog", 0.0}, {"sdlog", 1.0}}, {"meanlog", "sdlog"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Lognorm_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("lognorm", new_args_ptr, {"meanlog", "sdlog"}) { set_range({{"min", 0}}); }
    ~Lognorm_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Lognorm_Distribution : public Base_Distribution
{
public:
    Discrete_Lognorm_Distribution(dictargs new_args = {}) : Base_Distribution("lognorm_discr", {{"meanlog", 0.0}, {"sdlog", 1.0}}, {"meanlog", "sdlog"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Discrete_Lognorm_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("lognorm_discr", new_args_ptr, {"meanlog", "sdlog"}) { set_range({{"min", 0}}); }
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
    Gamma_Distribution(dictargs new_args = {}) : Base_Distribution("gamma", {{"shape", 1.0}, {"scale", 1.0}}, {"shape", "scale"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Gamma_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("gamma", new_args_ptr, {"shape", "scale"}) { set_range({{"min", 0}}); }
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
    Weibull_Distribution(dictargs new_args = {}) : Base_Distribution("weibull", {{"shape", 1.0}, {"scale", 1.0}}, {"shape", "scale"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Weibull_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("weibull", new_args_ptr, {"shape", "scale"}) { set_range({{"min", 0}}); }
    ~Weibull_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Weibull_Distribution : public Base_Distribution
{
public:
    Discrete_Weibull_Distribution(dictargs new_args = {}) : Base_Distribution("weibull_discr", {{"shape", 1.0}, {"scale", 1.0}}, {"shape", "scale"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Discrete_Weibull_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("weibull_discr", new_args_ptr, {"shape", "scale"}) { set_range({{"min", 0}}); }
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
    Rayleigh_Distribution(dictargs new_args = {}) : Base_Distribution("rayleigh", {{"scale", 1.0}}, {"scale"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Rayleigh_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("rayleigh", new_args_ptr, {"scale"}) { set_range({{"min", 0}}); }
    ~Rayleigh_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Rayleigh_Distribution : public Base_Distribution
{
public:
    Discrete_Rayleigh_Distribution(dictargs new_args = {}) : Base_Distribution("rayleigh_discr", {{"scale", 1.0}}, {"scale"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Discrete_Rayleigh_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("rayleigh_discr", new_args_ptr, {"scale"}) { set_range({{"min", 0}}); }
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
    Exponential_Distribution(dictargs new_args = {}) : Base_Distribution("exp", {{"rate", 1.0}}, {"rate"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Exponential_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("exp", new_args_ptr, {"rate"}) { set_range({{"min", 0}}); }
    ~Exponential_Distribution() {};

    virtual double pdf(double x) override;
    virtual double log_pdf(double x) override;
    virtual double logratio(double x, double y) override;
    virtual double rand(sys_rnd_gen &rng) override;
};

class Discrete_Exponential_Distribution : public Base_Distribution
{
public:
    Discrete_Exponential_Distribution(dictargs new_args = {}) : Base_Distribution("exp_discr", {{"rate", 1.0}}, {"rate"})
    {
        assign_args(new_args);
        set_range({{"min", 0}});
    }
    Discrete_Exponential_Distribution(dictargs_ptr new_args_ptr) : Base_Distribution("exp_discr", new_args_ptr, {"rate"}) { set_range({{"min", 0}}); }
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
template <typename dictargs_val_or_ptr>
Base_Distribution *distribution_ptr(const string &name = UNDEFINED, dictargs_val_or_ptr args = {})
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

#endif // DISTRIBUTIONS__H