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
/////////////////////////////  Uniform Distribution : methods  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Distr : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Distribution : public _Base_Distribution_
{
public:
    Uniform_Distribution(const dictargs &new_args = {}, const dictargs &range = {})
    {
        init("unif", {{"min", 0.0}, {"max", 1.0}}, {"min", "max"});
        assign_args(new_args);
        set_truncation(range);
    }
    Uniform_Distribution(const dictargs_ptr &new_args_ptr, const dictargs &range = {})
    {
        init("unif", new_args_ptr, {"min", "max"});
        set_truncation(range);
    }
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

    virtual double _cdf_(double x) const override
    {
        return x;
    }

    virtual double _logpdf_(double x) const override
    {
        UNUSED(x);
        return 0.0;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class Normal_Distribution : public _Base_Distribution_
{
public:
    Normal_Distribution(const dictargs &new_args = {}, const dictargs &range = {})
    {
        init("norm", {{"mean", 0.0}, {"sd", 1.0}}, {"mean", "sd"});
        assign_args(new_args);
        set_truncation(range);
    }
    Normal_Distribution(const dictargs_ptr &new_args_ptr, const dictargs &range = {})
    {
        init("norm", new_args_ptr, {"mean", "sd"});
        set_truncation(range);
    }
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
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            check_rand_attempts(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.39894228 * exp(-0.5 * pow((x - mean) / sd, 2)) / sd;
    }

    virtual double _cdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.5 + 0.5 * erf(0.707106781 * (x - mean) / sd);
    }

    virtual double _logpdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return -0.5 * pow((x - mean) / sd, 2) - 0.918938533 - log(sd);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Distribution : public _Base_Distribution_
{
public:
    Lognorm_Distribution(const dictargs &new_args = {}, const dictargs &range = {})
    {
        init("lnorm", {{"meanlog", 0.0}, {"sdlog", 1.0}}, {"meanlog", "sdlog"}, {{"min", 0}});
        assign_args(new_args);
        set_truncation(range);
    }
    Lognorm_Distribution(const dictargs_ptr &new_args_ptr, const dictargs &range = {})
    {
        init("lnorm", new_args_ptr, {"meanlog", "sdlog"}, {{"min", 0}});
        set_truncation(range);
    }
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
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            check_rand_attempts(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return exp(-0.5 * pow((log(x) - meanlog) / sdlog, 2)) * 0.39894228 / (x * sdlog);
    }

    virtual double _cdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return 0.5 + 0.5 * erf(0.707106781 * (log(x) - meanlog) / sdlog);
    }

    virtual double _logpdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return -0.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x * sdlog) - 0.918938533;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gamma_Distribution : public _Base_Distribution_
{
public:
    Gamma_Distribution(const dictargs &new_args = {}, const dictargs &range = {})
    {
        init("gamma", {{"shape", 1.0}, {"scale", 1.0}}, {"shape", "scale"}, {{"min", 0}});
        assign_args(new_args);
        set_truncation(range);
    }
    Gamma_Distribution(const dictargs_ptr &new_args_ptr, const dictargs &range = {})
    {
        init("gamma", new_args_ptr, {"shape", "scale"}, {{"min", 0}});
        set_truncation(range);
    }
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
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            check_rand_attempts(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
    }

    virtual double _cdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return auxgammap::gamma_cdf(x, shape, scale);
    }

    virtual double _logpdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Weibull_Distribution : public _Base_Distribution_
{
public:
    Weibull_Distribution(const dictargs &new_args = {}, const dictargs &range = {})
    {
        init("weibull", {{"shape", 1.0}, {"scale", 1.0}}, {"shape", "scale"}, {{"min", 0}});
        assign_args(new_args);
        set_truncation(range);
    }

    Weibull_Distribution(dictargs_ptr new_args_ptr, dictargs range = {})
    {
        init("weibull", new_args_ptr, {"shape", "scale"}, {{"min", 0}});
        set_truncation(range);
    }
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
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            check_rand_attempts(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
    }

    virtual double _cdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return 1. - exp(-pow(x / scale, shape));
    }

    virtual double _logpdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Rayleigh_Distr : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Rayleigh_Distribution : public _Base_Distribution_
{
public:
    Rayleigh_Distribution(const dictargs &new_args = {}, const dictargs &range = {})
    {
        init("rayleigh", {{"scale", 1.0}}, {"scale"}, {{"min", 0}});
        assign_args(new_args);
        set_truncation(range);
    }
    Rayleigh_Distribution(const dictargs_ptr &new_args_ptr, const dictargs &range = {})
    {
        init("rayleigh", new_args_ptr, {"scale"}, {{"min", 0}});
        set_truncation(range);
    }
    ~Rayleigh_Distribution() {};

    double logratio(double xold, double xnew) const override
    {
        double scale = get_arg(0);
        return log(xnew / xold) + 0.5 * pow(scale, -2) * (pow(xold, 2) - pow(xnew, 2));
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::uniform_real_distribution<double> unif_gen(0, 1);
        double xnew = get_arg(0) * sqrt(-2 * unif_gen(prng));
        int num_attempts = 0;
        while (!valid(xnew))
        {
            check_rand_attempts(num_attempts);
            xnew = get_arg(0) * sqrt(-2 * unif_gen(prng));
        }
        return xnew;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double scale = get_arg(0);
        return x * pow(scale, -2) * exp(-0.5 * pow(x / scale, 2));
    }

    virtual double _cdf_(double x) const override
    {
        double scale = get_arg(0);
        return exp(-0.5 * pow(x / scale, 2));
    }

    virtual double _logpdf_(double x) const override
    {
        double scale = get_arg(0);
        return log(x) - 2 * log(scale) - 0.5 * pow(x / scale, 2);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Exponential_Distr : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Exponential_Distribution : public _Base_Distribution_
{
public:
    Exponential_Distribution(dictargs new_args = {}, dictargs range = {})
    {
        init("exp", {{"rate", 1.0}}, {"rate"}, {{"min", 0}});
        assign_args(new_args);
        set_truncation(range);
    }
    Exponential_Distribution(const dictargs_ptr &new_args_ptr, const dictargs &range = {})
    {
        init("exp", new_args_ptr, {"rate"}, {{"min", 0}});
        set_truncation(range);
    }
    ~Exponential_Distribution() {};

    double logratio(double xold, double xnew) const override
    {
        double rate = get_arg(0);
        return rate * (xold - xnew);
    }

    double rand(sys_rnd_gen &prng) override
    {
        std::exponential_distribution<double> distr_gen(get_arg(0));
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            check_rand_attempts(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }

protected:
    virtual double _pdf_(double x) const override
    {
        double rate = get_arg(0);
        return rate * exp(-rate * x);
    }

    virtual double _cdf_(double x) const override
    {
        double rate = get_arg(0);
        return 1. - exp(-rate * x);
    }

    virtual double _logpdf_(double x) const override
    {
        double rate = get_arg(0);
        return log(rate) - rate * x;
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
std::shared_ptr<_Base_Distribution_> _distribution_ptr_(const std::string &name, const dictargs_val_or_ptr &args, const dictargs &range)
{
    std::shared_ptr<_Base_Distribution_> output = std::make_shared<_Base_Distribution_>();

    if (name == "unif")
    {
        output = std::make_shared<Uniform_Distribution>(args, range);
    }
    else if (name == "norm")
    {
        output = std::make_shared<Normal_Distribution>(args, range);
    }
    else if (name == "lnorm")
    {
        output = std::make_shared<Lognorm_Distribution>(args, range);
    }
    else if (name == "gamma")
    {
        output = std::make_shared<Gamma_Distribution>(args, range);
    }
    else if (name == "weibull")
    {
        output = std::make_shared<Weibull_Distribution>(args, range);
    }
    else if (name == "rayleigh")
    {
        output = std::make_shared<Rayleigh_Distribution>(args, range);
    }
    else if (name == "exp")
    {
        output = std::make_shared<Exponential_Distribution>(args, range);
    }
    else if (!tools::is_na(name))
    {
        tools::error("Unknown distribution '" + name + "'.");
    }

    return output;
}

/* Pointer to distribution with scalar arguments. */
std::shared_ptr<_Base_Distribution_> distribution_ptr(const std::string &name,
                                                      const dictargs &args,
                                                      const dictargs &range)
{
    return _distribution_ptr_(name, args, range);
}

/* Pointer to distribution with pointer arguments. */
std::shared_ptr<_Base_Distribution_> distribution_ptr(const std::string &name,
                                                      const dictargs_ptr &args,
                                                      const dictargs &range)
{
    return _distribution_ptr_(name, args, range);
}