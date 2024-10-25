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

/* A parent class to all distributions. */
class Base_Distribution : public Base_Function
{
public:
    using Base_Function::Base_Function;

    double pdf(double x) const { return _trunc_(_pdf_(x)); }
    double cdf(double x) const { return _trunc_(_shift_(_cdf_(x))); }
    double logpdf(double x) const { return _logtrunc_(_logpdf_(x)); }

    //----- Distribution specific methods -----//
    virtual double logratio(double xold, double xnew) const
    {
        UNUSED(xold);
        UNUSED(xnew);
        return 0.0;
    }
    virtual double rand(sys_rnd_gen &prng)
    {
        UNUSED(prng);
        return _undefined_call_error_("rand");
    }

protected:
    //----- Distribution specific methods -----//
    virtual double _pdf_(double x) const { return _undefined_call_error_("pdf", x); }
    virtual double _logpdf_(double x) const { return _undefined_call_error_("logpdf", x); }
    virtual double _cdf_(double x) const { return _undefined_call_error_("cdf", x); }

    //----- Adjust methods for truncation -----//
    void set_truncation(const dictargs &range);
    void check_rand_attempts(int &num_attempts) const;

private:
    //----- Internal error for calls to undefined -----//
    double _undefined_call_error_(const std::string &txt, double x = 0, double y = 0) const;


    //----- Functional methods for the truncation -----//
    double (Base_Distribution::*_trunc_method_)(double) const = &Base_Distribution::_trunc_dummy_;
    double _trunc_(double x) const { return (this->*_trunc_method_)(x); }
    double _trunc_dummy_(double x) const { return x; }
    double _trunc_min_(double x) const { return x / (1.0 - _cdf_(get_range_min())); }
    double _trunc_max_(double x) const { return x / (_cdf_(get_range_max())); }
    double _trunc_both_sides_(double x) const { return x / (_cdf_(get_range_max()) - _cdf_(get_range_max())); }

    //----- Functional methods for the (log) truncation -----//
    double (Base_Distribution::*_logtrunc_method_)(double) const = &Base_Distribution::_logtrunc_dummy_;
    double _logtrunc_(double x) const { return (this->*_logtrunc_method_)(x); }
    double _logtrunc_dummy_(double x) const { return x; }
    double _logtrunc_min_(double x) const { return x - log(1.0 - _cdf_(get_range_min())); }
    double _logtrunc_max_(double x) const { return x - log(_cdf_(get_range_max())); }
    double _logtrunc_both_sides_(double x) const { return x - log(_cdf_(get_range_max()) - _cdf_(get_range_min())); }

    //----- Functional methods to shift cumulative -----//
    double (Base_Distribution::*_shift_method_)(double) const = &Base_Distribution::_shift_dummy_;
    double _shift_(double x) const { return (this->*_shift_method_)(x); }
    double _shift_dummy_(double x) const { return x; }
    double _shift_min_(double x) const { return x - _cdf_(get_range_min()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Distribution : public Base_Distribution
{
public:
    Uniform_Distribution(
        dictargs new_args = {},
        dictargs range = {}) : Base_Distribution("unif",
                                                 {{"min", 0.0}, {"max", 1.0}},
                                                 {"min", "max"})
    {
        assign_args(new_args);
        set_truncation(range);
    }
    Uniform_Distribution(
        dictargs_ptr new_args_ptr,
        dictargs range = {}) : Base_Distribution("unif",
                                                 new_args_ptr,
                                                 {"min", "max"}) { set_truncation(range); }
    ~Uniform_Distribution() {}

    double logratio(double xold, double xnew) const override;
    double rand(sys_rnd_gen &prng) override;

protected:
    virtual double _pdf_(double x) const override;
    virtual double _logpdf_(double x) const override;
    virtual double _cdf_(double x) const override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Normal_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Distribution : public Base_Distribution
{
public:
    Normal_Distribution(
        dictargs new_args = {},
        dictargs range = {}) : Base_Distribution("norm",
                                                 {{"mean", 0.0}, {"sd", 1.0}},
                                                 {"mean", "sd"})
    {
        assign_args(new_args);
        set_truncation(range);
    }
    Normal_Distribution(
        dictargs_ptr new_args_ptr,
        dictargs range = {}) : Base_Distribution("norm",
                                                 new_args_ptr,
                                                 {"mean", "sd"}) { set_truncation(range); }
    ~Normal_Distribution() {};

    double logratio(double xold, double xnew) const override;
    double rand(sys_rnd_gen &prng) override;

protected:
    virtual double _pdf_(double x) const override;
    virtual double _logpdf_(double x) const override;
    virtual double _cdf_(double x) const override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Distribution : public Base_Distribution
{
public:
    Lognorm_Distribution(
        dictargs new_args = {},
        dictargs range = {}) : Base_Distribution("lnorm",
                                                 {{"meanlog", 0.0}, {"sdlog", 1.0}},
                                                 {"meanlog", "sdlog"},
                                                 {{"min", 0}})
    {
        assign_args(new_args);
        set_truncation(range);
    }
    Lognorm_Distribution(
        dictargs_ptr new_args_ptr,
        dictargs range = {}) : Base_Distribution("lnorm",
                                                 new_args_ptr,
                                                 {"meanlog", "sdlog"},
                                                 {{"min", 0}}) { set_truncation(range); }
    ~Lognorm_Distribution() {};

    double logratio(double xold, double xnew) const override;
    double rand(sys_rnd_gen &prng) override;

protected:
    virtual double _pdf_(double x) const override;
    virtual double _logpdf_(double x) const override;
    virtual double _cdf_(double x) const override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Gamma_Distr : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gamma_Distribution : public Base_Distribution
{
public:
    Gamma_Distribution(
        dictargs new_args = {},
        dictargs range = {}) : Base_Distribution("gamma",
                                                 {{"shape", 1.0}, {"scale", 1.0}},
                                                 {"shape", "scale"},
                                                 {{"min", 0}})
    {
        assign_args(new_args);
        set_truncation(range);
    }
    Gamma_Distribution(
        dictargs_ptr new_args_ptr,
        dictargs range = {}) : Base_Distribution("gamma",
                                                 new_args_ptr,
                                                 {"shape", "scale"},
                                                 {{"min", 0}}) { set_truncation(range); }
    ~Gamma_Distribution() {};

    double logratio(double xold, double xnew) const override;
    double rand(sys_rnd_gen &prng) override;

protected:
    virtual double _pdf_(double x) const override;
    virtual double _logpdf_(double x) const override;
    virtual double _cdf_(double x) const override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Weibull_Distr : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Weibull_Distribution : public Base_Distribution
{
public:
    Weibull_Distribution(
        dictargs new_args = {},
        dictargs range = {}) : Base_Distribution("weibull",
                                                 {{"shape", 1.0}, {"scale", 1.0}},
                                                 {"shape", "scale"},
                                                 {{"min", 0}})
    {
        assign_args(new_args);
        set_truncation(range);
    }

    Weibull_Distribution(
        dictargs_ptr new_args_ptr,
        dictargs range = {}) : Base_Distribution("weibull",
                                                 new_args_ptr,
                                                 {"shape", "scale"},
                                                 {{"min", 0}}) { set_truncation(range); }
    ~Weibull_Distribution() {};

    double logratio(double xold, double xnew) const override;
    double rand(sys_rnd_gen &prng) override;

protected:
    virtual double _pdf_(double x) const override;
    virtual double _logpdf_(double x) const override;
    virtual double _cdf_(double x) const override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Rayleigh_Distr : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Rayleigh_Distribution : public Base_Distribution
{
public:
    Rayleigh_Distribution(
        dictargs new_args = {},
        dictargs range = {}) : Base_Distribution("rayleigh",
                                                 {{"scale", 1.0}},
                                                 {"scale"},
                                                 {{"min", 0}})
    {
        assign_args(new_args);
        set_truncation(range);
    }
    Rayleigh_Distribution(
        dictargs_ptr new_args_ptr,
        dictargs range = {}) : Base_Distribution("rayleigh",
                                                 new_args_ptr,
                                                 {"scale"},
                                                 {{"min", 0}}) { set_truncation(range); }
    ~Rayleigh_Distribution() {};

    double logratio(double xold, double xnew) const override;
    double rand(sys_rnd_gen &prng) override;

protected:
    virtual double _pdf_(double x) const override;
    virtual double _logpdf_(double x) const override;
    virtual double _cdf_(double x) const override;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Exponential_Distr : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Exponential_Distribution : public Base_Distribution
{
public:
    Exponential_Distribution(
        dictargs new_args = {},
        dictargs range = {}) : Base_Distribution("exp",
                                                 {{"rate", 1.0}},
                                                 {"rate"},
                                                 {{"min", 0}})
    {
        assign_args(new_args);
        set_truncation(range);
    }
    Exponential_Distribution(
        dictargs_ptr new_args_ptr,
        dictargs range = {}) : Base_Distribution("exp",
                                                 new_args_ptr,
                                                 {"rate"},
                                                 {{"min", 0}}) { set_truncation(range); }
    ~Exponential_Distribution() {};

    double logratio(double xold, double xnew) const override;
    double rand(sys_rnd_gen &prng) override;

protected:
    virtual double _pdf_(double x) const override;
    virtual double _logpdf_(double x) const override;
    virtual double _cdf_(double x) const override;
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
Base_Distribution *distribution_ptr(const std::string &name = UNDEFINED, const dictargs_val_or_ptr &args = {}, const dictargs &range = {})
{
    Base_Distribution *output = new Base_Distribution();

    if (name == "unif")
    {
        output = new Uniform_Distribution(args, range);
    }
    else if (name == "norm")
    {
        output = new Normal_Distribution(args, range);
    }
    else if (name == "lnorm")
    {
        output = new Lognorm_Distribution(args, range);
    }
    else if (name == "gamma")
    {
        output = new Gamma_Distribution(args, range);
    }
    else if (name == "weibull")
    {
        output = new Weibull_Distribution(args, range);
    }
    else if (name == "rayleigh")
    {
        output = new Rayleigh_Distribution(args, range);
    }
    else if (name == "exp")
    {
        output = new Exponential_Distribution(args, range);
    }
    else if (!tools::is_na(name))
    {
        tools::error("Unknown distribution '" + name + "'.");
    }

    return output;
}

#endif // DISTRIBUTIONS__H