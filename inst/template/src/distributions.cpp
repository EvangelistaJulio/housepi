/*
 * distributions.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "varia.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Incomplete gamma Function  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME : this is all to avoid external libraries, but probably not worthy...
namespace auxgammap
{
    constexpr int GAMMAP_ITER = 10000;
    constexpr double GAMMAP_EPS = 1.e-12;
    constexpr double GAMMAP_FPMIN = 1e-30;

    /* Lower incomplete gamma function using series expansion. */
    double gamma_cdf_series(double scaled_x, double shape)
    {
        double sum = 1.0 / shape;
        double term = sum;
        for (int i = 1; (term > GAMMAP_EPS * sum) && (i < GAMMAP_ITER); ++i)
        {
            term *= scaled_x / (shape + i);
            sum += term;
        }
        return sum * exp(-scaled_x + shape * log(scaled_x) - lgamma(shape));
    }

    /* Upper incomplete gamma function using continued fraction. */
    double gamma_cdf_cf(double scaled_x, double shape)
    {
        double b = scaled_x + 1.0 - shape;
        double c = 1.0 / GAMMAP_FPMIN;
        double d = 1.0 / b;
        double h = d;
        double delta = 0.0;

        for (int i = 1; ((fabs(delta - 1.0) > GAMMAP_EPS) && (i < GAMMAP_ITER)); ++i)
        {
            double an = -i * (i - shape);
            b += 2.0;
            d = an * d + b;
            if (fabs(d) < GAMMAP_FPMIN)
            {
                d = GAMMAP_FPMIN;
            }
            c = b + an / c;
            if (fabs(c) < GAMMAP_FPMIN)
            {
                c = GAMMAP_FPMIN;
            }
            d = 1.0 / d;
            delta = c * d;
            h *= delta;
        }
        return 1.0 - exp(-scaled_x + shape * log(scaled_x) - lgamma(shape)) * h;
    }

    /* Function to compute the gamma cumulative distribution function. */
    double gamma_cdf(double x, double shape, double scale)
    {
        double scaled_x = x / scale;
        if (scaled_x < shape + 10.)
        {
            return gamma_cdf_series(scaled_x, shape);
        }
        else
        {
            return gamma_cdf_cf(scaled_x, shape);
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : methods  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Internal error for calls to undefined distribution. */
double Base_Distribution::_undefined_call_error_(const std::string &txt, double x, double y) const
{
    UNUSED(x);
    UNUSED(y);
    tools::error("Attempt to call method '" + txt + "' of undefined distribution.");
    return std::numeric_limits<double>::quiet_NaN();
}

/* Error for. */
void Base_Distribution::check_rand_attempts(int &num_attempts) const
{
    if (num_attempts > MAX_NUMBER_ATTEMPTS)
    {
        std::ostringstream oss;
        oss << "Exceeded attempts (N > " << MAX_NUMBER_ATTEMPTS << ")"
            << " in initialization of '" << get_name() << "' "
            << " within range {" << get_range_min() << ", " << get_range_max() << "}.";
        tools::error(oss.str());
    }
    ++num_attempts;
}

/* Truncate distribution at declaration.  */
void Base_Distribution::set_truncation(const dictargs &range)
{
    int trunc_cases = set_range(range);
    if (trunc_cases == 0) // no truncation
    {
        _trunc_method_ = &Base_Distribution::_trunc_dummy_;
        _logtrunc_method_ = &Base_Distribution::_logtrunc_dummy_;
        _shift_method_ = &Base_Distribution::_shift_dummy_;
    }
    else if (trunc_cases == 1) // lower range truncation
    {
        _trunc_method_ = &Base_Distribution::_trunc_min_;
        _logtrunc_method_ = &Base_Distribution::_logtrunc_min_;
        _shift_method_ = &Base_Distribution::_shift_min_;
    }
    else if (trunc_cases == 2) // upper range truncation
    {
        _trunc_method_ = &Base_Distribution::_trunc_max_;
        _logtrunc_method_ = &Base_Distribution::_logtrunc_max_;
        _shift_method_ = &Base_Distribution::_shift_dummy_;
    }
    else // both sides truncation
    {
        _trunc_method_ = &Base_Distribution::_trunc_both_sides_;
        _logtrunc_method_ = &Base_Distribution::_logtrunc_both_sides_;
        _shift_method_ = &Base_Distribution::_shift_dummy_;
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Uniform Distribution : methods  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

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

    double logratio(double xold, double xnew) const override
    {
        UNUSED(xold);
        UNUSED(xnew);
        return 0.0;
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::uniform_real_distribution<double> distr_gen(get_arg(0), get_arg(1));
        return distr_gen(prng);
    }

protected:
    virtual double _pdf_(double x) const override
    {
        UNUSED(x);
        return 1.0;
    }

    virtual double _logpdf_(double x) const override
    {
        UNUSED(x);
        return 0.0;
    }

    virtual double _cdf_(double x) const override
    {
        return x;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
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

    double logratio(double xold, double xnew) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.5 * (pow((xold - mean) / sd, 2) - pow((xnew - mean) / sd, 2));
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::normal_distribution<double> distr_gen(get_arg(0), get_arg(1));
        double new_value = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(new_value))
        {
            check_rand_attempts(num_attempts);
            new_value = distr_gen(prng);
        }
        return new_value;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.39894228 * exp(-0.5 * pow((x - mean) / sd, 2)) / sd;
    }

    virtual double _logpdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return -0.5 * pow((x - mean) / sd, 2) - 0.918938533 - log(sd);
    }

    virtual double _cdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.5 + 0.5 * erf(0.707106781 * (x - mean) / sd);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
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

    double logratio(double xold, double xnew) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return 0.5 * (pow((log(xold) - meanlog) / sdlog, 2) - pow((log(xnew) - meanlog) / sdlog, 2)) + log(xold / xnew);
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::lognormal_distribution<double> distr_gen(get_arg(0), get_arg(1));
        double new_value = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(new_value))
        {
            check_rand_attempts(num_attempts);
            new_value = distr_gen(prng);
        }
        return new_value;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return exp(-0.5 * pow((log(x) - meanlog) / sdlog, 2)) * 0.39894228 / (x * sdlog);
    }

    virtual double _logpdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return -0.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x * sdlog) - 0.918938533;
    }

    virtual double _cdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return 0.5 + 0.5 * erf(0.707106781 * (log(x) - meanlog) / sdlog);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
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

    double logratio(double xold, double xnew) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1) * (log(xnew / xold)) + (xold - xnew) / scale;
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::gamma_distribution<double> distr_gen(get_arg(0), get_arg(1));
        double new_value = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(new_value))
        {
            check_rand_attempts(num_attempts);
            new_value = distr_gen(prng);
        }
        return new_value;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
    }

    virtual double _logpdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
    }

    virtual double _cdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return auxgammap::gamma_cdf(x, shape, scale);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
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

    double logratio(double xold, double xnew) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1.) * (log(xnew / xold)) + pow(xold / scale, shape) - pow(xnew / scale, shape);
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::weibull_distribution<double> distr_gen(get_arg(0), get_arg(1));
        double new_value = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(new_value))
        {
            check_rand_attempts(num_attempts);
            new_value = distr_gen(prng);
        }
        return new_value;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
    }

    virtual double _logpdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
    }

    virtual double _cdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return 1. - exp(-pow(x / scale, shape));
    }
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

    double logratio(double xold, double xnew) const override
    {
        double scale = get_arg(0);
        return log(xnew / xold) + 0.5 * pow(scale, -2) * (pow(xold, 2) - pow(xnew, 2));
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::uniform_real_distribution<double> unif_gen(0, 1);
        double new_value = get_arg(0) * sqrt(-2 * unif_gen(prng));
        int num_attempts = 0;
        while (!valid(new_value))
        {
            check_rand_attempts(num_attempts);
            new_value = get_arg(0) * sqrt(-2 * unif_gen(prng));
        }
        return new_value;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double scale = get_arg(0);
        return x * pow(scale, -2) * exp(-0.5 * pow(x / scale, 2));
    }

    virtual double _logpdf_(double x) const override
    {
        double scale = get_arg(0);
        return log(x) - 2 * log(scale) - 0.5 * pow(x / scale, 2);
    }

    virtual double _cdf_(double x) const override
    {
        double scale = get_arg(0);
        return exp(-0.5 * pow(x / scale, 2));
    }
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

    double logratio(double xold, double xnew) const override
    {
        double rate = get_arg(0);
        return rate * (xold - xnew);
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::exponential_distribution<double> distr_gen(get_arg(0));
        double new_value = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(new_value))
        {
            check_rand_attempts(num_attempts);
            new_value = distr_gen(prng);
        }
        return new_value;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double rate = get_arg(0);
        return rate * exp(-rate * x);
    }

    virtual double _logpdf_(double x) const override
    {
        double rate = get_arg(0);
        return log(rate) - rate * x;
    }

    virtual double _cdf_(double x) const override
    {
        double rate = get_arg(0);
        return 1. - exp(-rate * x);
    }
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

/* Auxiliary hidden template to return the chosen distribution. */
template <typename dictargs_val_or_ptr>
Base_Distribution *_distribution_ptr_(const std::string &name, const dictargs_val_or_ptr &args, const dictargs &range)
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

/* Pointer to distribution with scalar arguments. */
Base_Distribution *distribution_ptr(const std::string &name,
                                    const dictargs &args,
                                    const dictargs &range)
{
    return _distribution_ptr_(name, args, range);
}

/* Pointer to distribution with pointer arguments. */
Base_Distribution *distribution_ptr(const std::string &name,
                                    const dictargs_ptr &args,
                                    const dictargs &range)
{
    return _distribution_ptr_(name, args, range);
}