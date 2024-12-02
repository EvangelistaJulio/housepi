/*
 * distribution.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "distribution.h"

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
////////////////////////////////  Base Distribution - wrapper  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::pdf(double x) const { return (this->*_func_pdf_)(x); }
double Base_Distribution::cdf(double x) const { return (this->*_func_cdf_)(x); }
double Base_Distribution::logpdf(double x) const { return (this->*_func_logpdf_)(x); }
double Base_Distribution::logratio(double xold, double xnew) const { return (this->*_func_logratio_)(xold, xnew); }
double Base_Distribution::rand(sys_rnd_gen &prng) { return (this->*_func_rand_)(prng); }

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution - virtual  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_pdf_(double x) const { return _undefined_call_error_("pdf", x); }
double Base_Distribution::_cdf_(double x) const { return _undefined_call_error_("cdf", x); }
double Base_Distribution::_logpdf_(double x) const { return _undefined_call_error_("logpdf", x); }
double Base_Distribution::_logratio_(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.0;
}
double Base_Distribution::_rand_(sys_rnd_gen &prng)
{
    UNUSED(prng);
    return _undefined_call_error_("rand");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : deaults  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_default_pdf_(double x) const { return (this->*_func_trunc_)(_pdf_(x)); }
double Base_Distribution::_default_cdf_(double x) const { return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(x))); }
double Base_Distribution::_default_logpdf_(double x) const { return (this->*_func_logtrunc_)(_logpdf_(x)); }
double Base_Distribution::_default_logratio_(double xold, double xnew) const { return _logratio_(xold, xnew); }
double Base_Distribution::_default_rand_(sys_rnd_gen &prng) { return _rand_(prng); }

/* Internal error for calls to undefined distribution. */
double Base_Distribution::_undefined_call_error_(const std::string &txt, double x, double y) const
{
    UNUSED(x);
    UNUSED(y);
    tools::error("Call to undefined method '" + txt + "' of distribution '" + get_name() + "'.");
    return std::numeric_limits<double>::quiet_NaN();
}

/* Check accept/reject within range is not exceeding maximum number of attempts. */
void Base_Distribution::_check_rand_attempts_(int &n_attempts) const
{
    if (n_attempts > MAX_NUMBER_ATTEMPTS)
    {
        std::ostringstream oss;
        oss << "Exceeded attempts (N > " << MAX_NUMBER_ATTEMPTS << ")"
            << " in random draw of " << get_name() << "(" << show_arg_values() << ")"
            << " with truncation {" << get_range_min() << ", " << get_range_max() << "}.";
        tools::error(oss.str());
    }
    ++n_attempts;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized center /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_discr_center_pdf_(double x) const
{
    double d = std::round(x);
    return (this->*_func_trunc_)(_cdf_(d + 0.5) - _cdf_(d - 0.5));
}

double Base_Distribution::_discr_center_cdf_(double x) const
{
    double d = std::round(x);
    return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(d + 0.5)));
}

double Base_Distribution::_discr_center_logpdf_(double x) const
{
    return log(_discr_center_pdf_(x));
}

double Base_Distribution::_discr_center_logratio_(double xold, double xnew) const
{
    double dold = std::round(xold), dnew = std::round(xnew);
    return log((_cdf_(dnew + 0.5) - _cdf_(dnew - 0.5)) / (_cdf_(dold + 0.5) - _cdf_(dold - 0.5)));
}

double Base_Distribution::_discr_center_rand_(sys_rnd_gen &prng)
{
    return std::round(_rand_(prng));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : discretized left /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_discr_left_pdf_(double x) const
{
    double d = std::floor(x);
    return (this->*_func_trunc_)(_cdf_(d + 1.0) - _cdf_(d));
}

double Base_Distribution::_discr_left_cdf_(double x) const
{
    double d = std::floor(x);
    return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(d + 1.0)));
}

double Base_Distribution::_discr_left_logpdf_(double x) const
{
    return log(_discr_left_pdf_(x));
}

double Base_Distribution::_discr_left_logratio_(double xold, double xnew) const
{
    double dold = std::floor(xold), dnew = std::floor(xnew);
    return log((_cdf_(dnew + 1.0) - _cdf_(dnew)) / (_cdf_(dold + 1.0) - _cdf_(dold)));
}

double Base_Distribution::_discr_left_rand_(sys_rnd_gen &prng)
{
    return std::floor(_rand_(prng));
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Base Distribution : truncation (functional)  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_trunc_dummy_(double x) const { return x; }
double Base_Distribution::_trunc_min_(double x) const { return x / (1.0 - _cdf_(get_range_min())); }
double Base_Distribution::_trunc_max_(double x) const { return x / (_cdf_(get_range_max())); }
double Base_Distribution::_trunc_both_sides_(double x) const { return x / (_cdf_(get_range_max()) - _cdf_(get_range_max())); }

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Base Distribution : logtruncation (functional)  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_logtrunc_dummy_(double x) const { return x; }
double Base_Distribution::_logtrunc_min_(double x) const { return x - log(1.0 - _cdf_(get_range_min())); }
double Base_Distribution::_logtrunc_max_(double x) const { return x - log(_cdf_(get_range_max())); }
double Base_Distribution::_logtrunc_both_sides_(double x) const { return x - log(_cdf_(get_range_max()) - _cdf_(get_range_min())); }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  Base Distribution : shift cumulative (functional)  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_shift_dummy_(double x) const { return x; }
double Base_Distribution::_shift_min_(double x) const { return x - _cdf_(get_range_min()); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : assign reshaping  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Truncate and/or discretize distribution at declaration.  */
void Base_Distribution::_assign_reshaping_(const dictargs &range, const std::string &discr)
{
    dictargs new_range = range;

    //===== DISCRETIZATION =====//
    // asjust range in case of discretization
    if (discr != UNDEFINED)
    {
        if (discr == "center")
        { // centered discretization
            _func_pdf_ = &Base_Distribution::_discr_center_pdf_;
            _func_cdf_ = &Base_Distribution::_discr_center_cdf_;
            _func_logpdf_ = &Base_Distribution::_discr_center_logpdf_;
            _func_logratio_ = &Base_Distribution::_discr_center_logratio_;
            _func_rand_ = &Base_Distribution::_discr_center_rand_;

            if (range.find("min") != range.end())
            {
                new_range["min"] = std::floor(range.at("min")) + 0.5;
            }
            if (range.find("max") != range.end())
            {
                new_range["max"] = std::ceil(range.at("max")) - 0.5;
            }
        }
        else if (discr == "left")
        { // left adjusted discretization
            _func_pdf_ = &Base_Distribution::_discr_left_pdf_;
            _func_cdf_ = &Base_Distribution::_discr_left_cdf_;
            _func_logpdf_ = &Base_Distribution::_discr_left_logpdf_;
            _func_logratio_ = &Base_Distribution::_discr_left_logratio_;
            _func_rand_ = &Base_Distribution::_discr_left_rand_;

            if (range.find("min") != range.end())
            {
                new_range["min"] = std::floor(range.at("min"));
            }
            if (range.find("max") != range.end())
            {
                new_range["max"] = std::ceil(range.at("max"));
            }
        }
        else
        {
            tools::error(
                "Unrecognized discretization method '" + discr + "'.");
        }
    }

    //===== TRUNCATION =====/
    int trunc_cases = set_range(new_range);

    //----- Fucntional methods -----//
    if (trunc_cases == 0) // no truncation
    {
        _func_trunc_ = &Base_Distribution::_trunc_dummy_;
        _func_logtrunc_ = &Base_Distribution::_logtrunc_dummy_;
        _func_shift_ = &Base_Distribution::_shift_dummy_;
    }
    else if (trunc_cases == 1) // lower range truncation
    {
        _func_trunc_ = &Base_Distribution::_trunc_min_;
        _func_logtrunc_ = &Base_Distribution::_logtrunc_min_;
        _func_shift_ = &Base_Distribution::_shift_min_;
    }
    else if (trunc_cases == 2) // upper range truncation
    {
        _func_trunc_ = &Base_Distribution::_trunc_max_;
        _func_logtrunc_ = &Base_Distribution::_logtrunc_max_;
        _func_shift_ = &Base_Distribution::_shift_dummy_;
    }
    else // both sides truncation
    {
        _func_trunc_ = &Base_Distribution::_trunc_both_sides_;
        _func_logtrunc_ = &Base_Distribution::_logtrunc_both_sides_;
        _func_shift_ = &Base_Distribution::_shift_dummy_;
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

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
        int n_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(n_attempts);
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
        int n_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(n_attempts);
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
        int n_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(n_attempts);
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
        int n_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(n_attempts);
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
        int n_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(n_attempts);
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
        int n_attempts = 0;
        while (!valid(xnew))
        {
            _check_rand_attempts_(n_attempts);
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
std::shared_ptr<Base_Distribution> _distribution_switchboard_(const std::string &name, const dictargs_val_or_ptr &args, const dictargs &range)
{
    std::shared_ptr<Base_Distribution> output = std::make_shared<Base_Distribution>();
    const std::vector<std::string> options = {
        "unif",
        "norm",
        "lnorm",
        "gamma",
        "weibull",
        "rayleigh",
        "exp"};

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
        std::string err_message = "Unknown distribution '" + name + "'.";
        err_message += tools::find_closest_option(name, options);
        throw std::invalid_argument(err_message);
    }

    return output;
}

/* Pointer to distribution with scalar arguments. */
std::shared_ptr<Base_Distribution> distribution_switchboard(const std::string &name,
                                                            const dictargs &args,
                                                            const dictargs &range)
{
    return _distribution_switchboard_(name, args, range);
}

/* Pointer to distribution with pointer arguments. */
std::shared_ptr<Base_Distribution> distribution_switchboard(const std::string &name,
                                                            const dictargs_ptr &args,
                                                            const dictargs &range)
{
    return _distribution_switchboard_(name, args, range);
}
