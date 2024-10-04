/*
 * distributions.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef DISTRIBUTIONS__H
#define DISTRIBUTIONS__H

#include <random>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

#include "utils.h"

// FIXME : this is all to avoid external libraries, but probably not worthy...
namespace auxgammap
{
    constexpr int GAMMAP_ITER = 10000;
    constexpr double GAMMAP_EPS = 1.e-12;
    constexpr double GAMMAP_FPMIN = 1e-30;

    double gamma_cdf_series(double scaled_x, double shape);
    double gamma_cdf_cf(double scaled_x, double shape);
    double gamma_cdf(double x, double shape, double scale);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Base_Distribution : public Base_Function
{
public:
    using Base_Function::Base_Function;

    virtual double pdf(double x) const
    {
        UNUSED(x);
        _undefined_call_error_();
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double cdf(double x) const
    {
        UNUSED(x);
        _undefined_call_error_();
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double logpdf(double x) const
    {
        UNUSED(x);
        _undefined_call_error_();
        return std::numeric_limits<double>::quiet_NaN();
    }
    virtual double logratio(double x, double y) const
    {
        UNUSED(x);
        UNUSED(y);
        return 0.0;
    }
    virtual double rand(sys_rnd_gen &rng)
    {
        UNUSED(rng);
        _undefined_call_error_();
        return std::numeric_limits<double>::quiet_NaN();
    }

private:
    void _undefined_call_error_() const
    {
        basics::raise_error("Attempt to call a method of an undefined distribution.");
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Distribution : public Base_Distribution
{
public:
    Uniform_Distribution(
        dictargs new_args = {}) : Base_Distribution("unif",
                                                    {{"min", 0.0}, {"max", 1.0}},
                                                    {"min", "max"}) { assign_args(new_args); }
    Uniform_Distribution(
        dictargs_ptr new_args_ptr) : Base_Distribution("unif",
                                                       new_args_ptr,
                                                       {"min", "max"}) {}
    ~Uniform_Distribution() {}

    double pdf(double x) const override;
    double cdf(double x) const override;
    double logpdf(double x) const override;
    double logratio(double x, double y) const override;
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Normal_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Distribution : public Base_Distribution
{
public:
    Normal_Distribution(
        dictargs new_args = {}) : Base_Distribution("norm",
                                                    {{"mean", 0.0}, {"sd", 1.0}},
                                                    {"mean", "sd"}) { assign_args(new_args); }
    Normal_Distribution(
        dictargs_ptr new_args_ptr) : Base_Distribution("norm",
                                                       new_args_ptr,
                                                       {"mean", "sd"}) {}
    ~Normal_Distribution() {};

    double pdf(double x) const override;
    double cdf(double x) const override;
    double logpdf(double x) const override;
    double logratio(double x, double y) const override;
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Distribution : public Base_Distribution
{
public:
    Lognorm_Distribution(
        dictargs new_args = {}) : Base_Distribution("lnorm",
                                                    {{"meanlog", 0.0}, {"sdlog", 1.0}},
                                                    {"meanlog", "sdlog"},
                                                    {{"min", 0}}) { assign_args(new_args); }
    Lognorm_Distribution(
        dictargs_ptr new_args_ptr) : Base_Distribution("lnorm",
                                                       new_args_ptr,
                                                       {"meanlog", "sdlog"},
                                                       {{"min", 0}}) {}
    ~Lognorm_Distribution() {};

    double pdf(double x) const override;
    double cdf(double x) const override;
    double logpdf(double x) const override;
    double logratio(double x, double y) const override;
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Gamma_Distr : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gamma_Distribution : public Base_Distribution
{
public:
    Gamma_Distribution(
        dictargs new_args = {}) : Base_Distribution("gamma",
                                                    {{"shape", 1.0}, {"scale", 1.0}},
                                                    {"shape", "scale"},
                                                    {{"min", 0}}) { assign_args(new_args); }
    Gamma_Distribution(
        dictargs_ptr new_args_ptr) : Base_Distribution("gamma",
                                                       new_args_ptr,
                                                       {"shape", "scale"},
                                                       {{"min", 0}}) {}
    ~Gamma_Distribution() {};

    double pdf(double x) const override;
    double cdf(double x) const override;
    double logpdf(double x) const override;
    double logratio(double x, double y) const override;
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Weibull_Distr : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Weibull_Distribution : public Base_Distribution
{
public:
    Weibull_Distribution(
        dictargs new_args = {}) : Base_Distribution("weibull",
                                                    {{"shape", 1.0}, {"scale", 1.0}},
                                                    {"shape", "scale"},
                                                    {{"min", 0}}) { assign_args(new_args); }

    Weibull_Distribution(
        dictargs_ptr new_args_ptr) : Base_Distribution("weibull",
                                                       new_args_ptr,
                                                       {"shape", "scale"},
                                                       {{"min", 0}}) {}
    ~Weibull_Distribution() {};

    double pdf(double x) const override;
    double cdf(double x) const override;
    double logpdf(double x) const override;
    double logratio(double x, double y) const override;
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Rayleigh_Distr : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Rayleigh_Distribution : public Base_Distribution
{
public:
    Rayleigh_Distribution(
        dictargs new_args = {}) : Base_Distribution("rayleigh",
                                                    {{"scale", 1.0}},
                                                    {"scale"},
                                                    {{"min", 0}}) { assign_args(new_args); }
    Rayleigh_Distribution(
        dictargs_ptr new_args_ptr) : Base_Distribution("rayleigh",
                                                       new_args_ptr,
                                                       {"scale"},
                                                       {{"min", 0}}) {}
    ~Rayleigh_Distribution() {};

    double pdf(double x) const override;
    double cdf(double x) const override;
    double logpdf(double x) const override;
    double logratio(double x, double y) const override;
    double rand(sys_rnd_gen &rng) override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Exponential_Distr : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Exponential_Distribution : public Base_Distribution
{
public:
    Exponential_Distribution(
        dictargs new_args = {}) : Base_Distribution("exp",
                                                    {{"rate", 1.0}},
                                                    {"rate"},
                                                    {{"min", 0}}) { assign_args(new_args); }
    Exponential_Distribution(
        dictargs_ptr new_args_ptr) : Base_Distribution("exp",
                                                       new_args_ptr,
                                                       {"rate"},
                                                       {{"min", 0}}) {}
    ~Exponential_Distribution() {};

    double pdf(double x) const override;
    double cdf(double x) const override;
    double logpdf(double x) const override;
    double logratio(double x, double y) const override;
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

    if (name == "unif")
    {
        output = new Uniform_Distribution(args);
    }
    else if (name == "norm")
    {
        output = new Normal_Distribution(args);
    }
    else if (name == "lnorm")
    {
        output = new Lognorm_Distribution(args);
    }
    else if (name == "gamma")
    {
        output = new Gamma_Distribution(args);
    }
    else if (name == "weibull")
    {
        output = new Weibull_Distribution(args);
    }
    else if (name == "rayleigh")
    {
        output = new Rayleigh_Distribution(args);
    }
    else if (name == "exp")
    {
        output = new Exponential_Distribution(args);
    }
    else if (!basics::is_na(name))
    {
        basics::raise_error("Unknown distribution '" + name + "'.");
    }

    return output;
}

#endif // DISTRIBUTIONS__H