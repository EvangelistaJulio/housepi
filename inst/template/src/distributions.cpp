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

class Uniform_Distribution : public Base_Distribution
{
private:
    std::uniform_real_distribution<double> distr_gen;

public:
    Uniform_Distribution(const dictargs &new_args = {},
                         const dictargs &range = {},
                         bool discr = false) : distr_gen(0.0, 1.0)
    {
        init("unif", {{"min", 0.0}, {"max", 1.0}}, {"min", "max"});
        assign_args(new_args);
        _assign_reshaping_(range, ((discr) ? "center" : UNDEFINED));
    }

    Uniform_Distribution(const dictargs_ptr &new_args_ptr,
                         const dictargs &range = {},
                         bool discr = false) : distr_gen(0.0, 1.0)
    {
        init("unif", new_args_ptr, {"min", "max"});
        _assign_reshaping_(range, ((discr) ? "center" : UNDEFINED));
    }

    ~Uniform_Distribution() {}

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

    virtual double _logratio_(double xold, double xnew) const override
    {
        UNUSED(xold);
        UNUSED(xnew);
        return 0.0;
    }

    virtual double _rand_(sys_rnd_gen &prng) override
    {
        distr_gen.param(std::uniform_real_distribution<double>::param_type(get_arg(0), get_arg(1)));
        return distr_gen(prng);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
class Normal_Distribution : public Base_Distribution
{
private:
    std::normal_distribution<double> distr_gen;

public:
    Normal_Distribution(const dictargs &new_args = {},
                        const dictargs &range = {},
                        bool discr = false) : distr_gen(0.0, 1.0)
    {
        init("norm", {{"mean", 0.0}, {"sd", 1.0}}, {"mean", "sd"});
        assign_args(new_args);
        _assign_reshaping_(range, ((discr) ? "center" : UNDEFINED));
    }

    Normal_Distribution(const dictargs_ptr &new_args_ptr,
                        const dictargs &range = {},
                        bool discr = false) : distr_gen(0.0, 1.0)
    {
        init("norm", new_args_ptr, {"mean", "sd"});
        _assign_reshaping_(range, ((discr) ? "center" : UNDEFINED));
    }

    ~Normal_Distribution() {};

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

    virtual double _logratio_(double xold, double xnew) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.5 * (pow((xold - mean) / sd, 2) - pow((xnew - mean) / sd, 2));
    }

    virtual double _rand_(sys_rnd_gen &prng) override
    {
        distr_gen.param(std::normal_distribution<double>::param_type(get_arg(0), get_arg(1)));
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Distribution : public Base_Distribution
{
private:
    std::lognormal_distribution<double> distr_gen;

public:
    Lognorm_Distribution(const dictargs &new_args = {},
                         const dictargs &range = {},
                         bool discr = false) : distr_gen(0.0, 1.0)
    {
        init("lnorm", {{"meanlog", 0.0}, {"sdlog", 1.0}}, {"meanlog", "sdlog"}, {{"min", 0}});
        assign_args(new_args);
        _assign_reshaping_(range, ((discr) ? "center" : UNDEFINED));
    }

    Lognorm_Distribution(const dictargs_ptr &new_args_ptr,
                         const dictargs &range = {},
                         bool discr = false) : distr_gen(0.0, 1.0)
    {
        init("lnorm", new_args_ptr, {"meanlog", "sdlog"}, {{"min", 0}});
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }

    ~Lognorm_Distribution() {};

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

    virtual double _logratio_(double xold, double xnew) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return 0.5 * (pow((log(xold) - meanlog) / sdlog, 2) - pow((log(xnew) - meanlog) / sdlog, 2)) + log(xold / xnew);
    }

    virtual double _rand_(sys_rnd_gen &prng) override
    {
        distr_gen.param(std::lognormal_distribution<double>::param_type(get_arg(0), get_arg(1)));
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Gamma_Distribution : public Base_Distribution
{
private:
    std::gamma_distribution<double> distr_gen;

public:
    Gamma_Distribution(const dictargs &new_args = {},
                       const dictargs &range = {},
                       bool discr = false) : distr_gen(1.0, 1.0)
    {
        init("gamma", {{"shape", 1.0}, {"scale", 1.0}}, {"shape", "scale"}, {{"min", 0}});
        assign_args(new_args);
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }

    Gamma_Distribution(const dictargs_ptr &new_args_ptr,
                       const dictargs &range = {},
                       bool discr = false) : distr_gen(1.0, 1.0)
    {
        init("gamma", new_args_ptr, {"shape", "scale"}, {{"min", 0}});
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }

    ~Gamma_Distribution() {};

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

    virtual double _logratio_(double xold, double xnew) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1) * (log(xnew / xold)) + (xold - xnew) / scale;
    }

    virtual double _rand_(sys_rnd_gen &prng) override
    {
        distr_gen.param(std::gamma_distribution<double>::param_type(get_arg(0), get_arg(1)));
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Weibull_Distribution : public Base_Distribution
{
private:
    std::weibull_distribution<double> distr_gen;

public:
    Weibull_Distribution(const dictargs &new_args = {},
                         const dictargs &range = {},
                         bool discr = false) : distr_gen(1.0, 1.0)
    {
        init("weibull", {{"shape", 1.0}, {"scale", 1.0}}, {"shape", "scale"}, {{"min", 0}});
        assign_args(new_args);
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }

    Weibull_Distribution(const dictargs_ptr new_args_ptr,
                         const dictargs &range = {},
                         bool discr = false) : distr_gen(1.0, 1.0)
    {
        init("weibull", new_args_ptr, {"shape", "scale"}, {{"min", 0}});
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }
    ~Weibull_Distribution() {};

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

    virtual double _logratio_(double xold, double xnew) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1.) * (log(xnew / xold)) + pow(xold / scale, shape) - pow(xnew / scale, shape);
    }

    virtual double _rand_(sys_rnd_gen &prng) override
    {
        distr_gen.param(std::weibull_distribution<double>::param_type(get_arg(0), get_arg(1)));
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Rayleigh_Distr : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Rayleigh_Distribution : public Base_Distribution
{
private:
    std::uniform_real_distribution<double> unif_gen;

public:
    Rayleigh_Distribution(const dictargs &new_args = {},
                          const dictargs &range = {},
                          bool discr = false) : unif_gen(0.0, 1.0)
    {
        init("rayleigh", {{"scale", 1.0}}, {"scale"}, {{"min", 0}});
        assign_args(new_args);
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }
    Rayleigh_Distribution(const dictargs_ptr &new_args_ptr,
                          const dictargs &range = {},
                          bool discr = false) : unif_gen(0.0, 1.0)
    {
        init("rayleigh", new_args_ptr, {"scale"}, {{"min", 0}});
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }
    ~Rayleigh_Distribution() {};

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

    virtual double _logratio_(double xold, double xnew) const override
    {
        double scale = get_arg(0);
        return log(xnew / xold) + 0.5 * pow(scale, -2) * (pow(xold, 2) - pow(xnew, 2));
    }

    virtual double _rand_(sys_rnd_gen &prng) override
    {
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
    }

private:
    double distr_gen(sys_rnd_gen &prng)
    {
        return get_arg(0) * sqrt(-2 * unif_gen(prng));
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Exponential_Distr : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Exponential_Distribution : public Base_Distribution
{
private:
    std::exponential_distribution<double> distr_gen;

public:
    Exponential_Distribution(const dictargs &new_args = {},
                             const dictargs &range = {},
                             bool discr = false) : distr_gen(1.0)
    {
        init("exp", {{"rate", 1.0}}, {"rate"}, {{"min", 0}});
        assign_args(new_args);
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }
    Exponential_Distribution(const dictargs_ptr &new_args_ptr,
                             const dictargs &range = {},
                             bool discr = false) : distr_gen(1.0)
    {
        init("exp", new_args_ptr, {"rate"}, {{"min", 0}});
        _assign_reshaping_(range, ((discr) ? "left" : UNDEFINED));
    }
    ~Exponential_Distribution() {};

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

    virtual double _logratio_(double xold, double xnew) const override
    {
        double rate = get_arg(0);
        return rate * (xold - xnew);
    }

    virtual double _rand_(sys_rnd_gen &prng) override
    {
        distr_gen.param(std::exponential_distribution<double>::param_type(get_arg(0)));
        double xnew = distr_gen(prng);
        int num_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(num_attempts);
            xnew = distr_gen(prng);
        }
        return xnew;
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
std::shared_ptr<Base_Distribution> _distribution_ptr_(const std::string &name, const dictargs_val_or_ptr &args, const dictargs &range)
{
    std::shared_ptr<Base_Distribution> output = std::make_shared<Base_Distribution>();

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
std::shared_ptr<Base_Distribution> distribution_ptr(const std::string &name,
                                                      const dictargs &args,
                                                      const dictargs &range)
{
    return _distribution_ptr_(name, args, range);
}

/* Pointer to distribution with pointer arguments. */
std::shared_ptr<Base_Distribution> distribution_ptr(const std::string &name,
                                                      const dictargs_ptr &args,
                                                      const dictargs &range)
{
    return _distribution_ptr_(name, args, range);
}