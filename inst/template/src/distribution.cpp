/*
 * distribution.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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
double Base_Distribution::rand(rnd_gen &prng)
{
    _update_args_();
    return (this->*_func_rand_)(prng);
}

void Base_Distribution::assign_reshaping(bool discrete,
                                         double lower,
                                         double upper)
{
    _assign_reshaping_((discrete) ? "center" : NULL_STRING,
                       lower, upper);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution - virtual  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_pdf_(double x) const { return _undefined_call_("pdf", x); }
double Base_Distribution::_cdf_(double x) const { return _undefined_call_("cdf", x); }
double Base_Distribution::_logpdf_(double x) const { return _undefined_call_("logpdf", x); }
double Base_Distribution::_logratio_(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.0;
}
void Base_Distribution::_update_args_() {};
double Base_Distribution::_rand_(rnd_gen &prng)
{
    UNUSED(prng);
    return _undefined_call_("rand");
}

/* Internal function for calls to undefined distribution. */
double Base_Distribution::_undefined_call_(const std::string &txt, double x, double y) const
{
    UNUSED(x);
    UNUSED(y);
    throw std::invalid_argument(
        "Call to undefined method \"" + txt + "\" of distribution \"" + get_name() + "\".");

    return std::numeric_limits<double>::quiet_NaN();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : deaults  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Distribution::_default_pdf_(double x) const { return (this->*_func_trunc_)(_pdf_(x)); }
double Base_Distribution::_default_cdf_(double x) const { return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(x))); }
double Base_Distribution::_default_logpdf_(double x) const { return (this->*_func_logtrunc_)(_logpdf_(x)); }
double Base_Distribution::_default_logratio_(double xold, double xnew) const { return _logratio_(xold, xnew); }

/* Check accept/reject within range is not exceeding maximum number of attempts. */
double Base_Distribution::_default_rand_(rnd_gen &prng)
{
    double xnew = _rand_(prng);

    unsigned n_attempts = 0;
    while (!valid(xnew))
    {
        if (n_attempts > MAX_NUMBER_ATTEMPTS)
        {
            std::ostringstream oss;
            oss << "Exceeded attempts (N > " << MAX_NUMBER_ATTEMPTS << ")"
                << " in random draw of " << get_name() << "(" << show_arg_values() << ")"
                << " with truncation {" << get_range_min() << ", " << get_range_max() << "}.";
            std::invalid_argument(oss.str());
        }
        xnew = _rand_(prng);
        ++n_attempts;
    }

    return xnew;
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

double Base_Distribution::_discr_center_rand_(rnd_gen &prng)
{
    // note : lround ro avoid casting like round(-0.1) -> -0.0
    return static_cast<double>(std::lround(_default_rand_(prng)));
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

double Base_Distribution::_discr_left_rand_(rnd_gen &prng)
{
    return std::floor(_default_rand_(prng));
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
void Base_Distribution::_assign_reshaping_(const std::string &discrete,
                                           double lower,
                                           double upper)
{
    //===== DISCRETIZATION =====//
    // asjust range in case of discretization
    if (!tools::is_na(discrete))
    {
        if (discrete == "center")
        { // centered discretization
            _func_pdf_ = &Base_Distribution::_discr_center_pdf_;
            _func_cdf_ = &Base_Distribution::_discr_center_cdf_;
            _func_logpdf_ = &Base_Distribution::_discr_center_logpdf_;
            _func_logratio_ = &Base_Distribution::_discr_center_logratio_;
            _func_rand_ = &Base_Distribution::_discr_center_rand_;

            if (!tools::is_na(lower))
                lower = std::floor(lower) - 0.5;

            if (!tools::is_na(upper))
                upper = std::ceil(upper) + 0.5;
        }
        else if (discrete == "left")
        { // left adjusted discretization
            _func_pdf_ = &Base_Distribution::_discr_left_pdf_;
            _func_cdf_ = &Base_Distribution::_discr_left_cdf_;
            _func_logpdf_ = &Base_Distribution::_discr_left_logpdf_;
            _func_logratio_ = &Base_Distribution::_discr_left_logratio_;
            _func_rand_ = &Base_Distribution::_discr_left_rand_;

            if (!tools::is_na(lower))
                lower = std::floor(lower);

            if (!tools::is_na(upper))
                upper = std::ceil(upper);
        }
        else
        {
            throw std::invalid_argument(
                "Unrecognized discretization method \"" + discrete + "\".");
        }
    }

    //===== TRUNCATION =====/
    ushort trunc_cases = set_range(lower, upper);

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
    Uniform_Distribution() : distr_gen(0.0, 1.0)
    {
        _init_defaults_("unif", {"min", "max"}, {0.0, 1.0},
                        NULL_DOUBLE, NULL_DOUBLE, false);
    }

    ~Uniform_Distribution() {}

protected:
    double _pdf_(double x) const override
    {
        UNUSED(x);
        return 1.0;
    }

    double _cdf_(double x) const override
    {
        return x;
    }

    double _logpdf_(double x) const override
    {
        UNUSED(x);
        return 0.0;
    }

    double _logratio_(double xold, double xnew) const override
    {
        UNUSED(xold);
        UNUSED(xnew);
        return 0.0;
    }

    void _update_args_() override
    {
        distr_gen.param(std::uniform_real_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    double _rand_(rnd_gen &prng) override
    {
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
    Normal_Distribution() : distr_gen(0.0, 1.0)
    {
        _init_defaults_("norm", {"mean", "sd"}, {0.0, 1.0},
                        NULL_DOUBLE, NULL_DOUBLE, false);
    }

    ~Normal_Distribution() {};

protected:
    double _pdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.39894228 * exp(-0.5 * pow((x - mean) / sd, 2)) / sd;
    }

    double _cdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.5 + 0.5 * erf(0.707106781 * (x - mean) / sd);
    }

    double _logpdf_(double x) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return -0.5 * pow((x - mean) / sd, 2) - 0.918938533 - log(sd);
    }

    double _logratio_(double xold, double xnew) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        return 0.5 * (pow((xold - mean) / sd, 2) - pow((xnew - mean) / sd, 2));
    }

    void _update_args_() override
    {
        distr_gen.param(std::normal_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    double _rand_(rnd_gen &prng) override
    {
        return distr_gen(prng);
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
    Lognorm_Distribution() : distr_gen(0.0, 1.0)
    {
        _init_defaults_("lnorm",
                        {"meanlog", "sdlog"}, {0.0, 1.0},
                        0.0, NULL_DOUBLE, true);
    }
    ~Lognorm_Distribution() {};

    void assign_reshaping(bool discrete, double lower, double upper) override
    {
        if ((!tools::is_na(lower)) && (lower <= 0))
            throw std::invalid_argument(
                "Invalid min value " + std::to_string(lower) + " for distribution \"" + get_name() + "\".");
        _assign_reshaping_(((discrete) ? "left" : NULL_STRING), lower, upper);
    }

protected:
    double _pdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return exp(-0.5 * pow((log(x) - meanlog) / sdlog, 2)) * 0.39894228 / (x * sdlog);
    }

    double _cdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return 0.5 + 0.5 * erf(0.707106781 * (log(x) - meanlog) / sdlog);
    }

    double _logpdf_(double x) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return -0.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x * sdlog) - 0.918938533;
    }

    double _logratio_(double xold, double xnew) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        return 0.5 * (pow((log(xold) - meanlog) / sdlog, 2) - pow((log(xnew) - meanlog) / sdlog, 2)) + log(xold / xnew);
    }

    void _update_args_() override
    {
        distr_gen.param(std::lognormal_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    double _rand_(rnd_gen &prng) override
    {
        return distr_gen(prng);
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
    Gamma_Distribution() : distr_gen(1.0, 1.0)
    {
        _init_defaults_("gamma",
                        {"shape", "scale"}, {1.0, 1.0},
                        0.0, NULL_DOUBLE, true);
    }
    ~Gamma_Distribution() {};

    void assign_reshaping(bool discrete, double lower, double upper) override
    {
        if ((!tools::is_na(lower)) && (lower <= 0))
            throw std::invalid_argument(
                "Invalid min value " + std::to_string(lower) + " for distribution \"" + get_name() + "\".");
        _assign_reshaping_(((discrete) ? "left" : NULL_STRING), lower, upper);
    }

protected:
    double _pdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
    }

    double _cdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return auxgammap::gamma_cdf(x, shape, scale);
    }

    double _logpdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
    }

    double _logratio_(double xold, double xnew) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1) * (log(xnew / xold)) + (xold - xnew) / scale;
    }

    void _update_args_() override
    {
        distr_gen.param(std::gamma_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    double _rand_(rnd_gen &prng) override
    {
        return distr_gen(prng);
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
    Weibull_Distribution() : distr_gen(1.0, 1.0)
    {
        _init_defaults_("weibull",
                        {"shape", "scale"}, {1.0, 1.0},
                        0.0, NULL_DOUBLE, true);
    }
    ~Weibull_Distribution() {};

    void assign_reshaping(bool discrete, double lower, double upper) override
    {
        if ((!tools::is_na(lower)) && (lower <= 0))
            throw std::invalid_argument(
                "Invalid min value " + std::to_string(lower) + " for distribution \"" + get_name() + "\".");
        _assign_reshaping_(((discrete) ? "left" : NULL_STRING), lower, upper);
    }

protected:
    double _pdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
    }

    double _cdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return 1. - exp(-pow(x / scale, shape));
    }

    double _logpdf_(double x) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
    }

    double _logratio_(double xold, double xnew) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        return (shape - 1.) * (log(xnew / xold)) + pow(xold / scale, shape) - pow(xnew / scale, shape);
    }

    void _update_args_() override
    {
        distr_gen.param(std::weibull_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    double _rand_(rnd_gen &prng) override
    {
        return distr_gen(prng);
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
    Rayleigh_Distribution() : unif_gen(0.0, 1.0)
    {
        _init_defaults_("rayleigh",
                        {"scale"}, {1.0},
                        0.0, NULL_DOUBLE, true);
    }
    ~Rayleigh_Distribution() {};

    void assign_reshaping(bool discrete, double lower, double upper) override
    {
        if ((!tools::is_na(lower)) && (lower <= 0))
            throw std::invalid_argument(
                "Invalid min value " + std::to_string(lower) + " for distribution \"" + get_name() + "\".");
        _assign_reshaping_(((discrete) ? "left" : NULL_STRING), lower, upper);
    }

protected:
    double _pdf_(double x) const override
    {
        double scale = get_arg(0);
        return x * pow(scale, -2) * exp(-0.5 * pow(x / scale, 2));
    }

    double _cdf_(double x) const override
    {
        double scale = get_arg(0);
        return 1. - exp(-0.5 * pow(x / scale, 2));
    }

    double _logpdf_(double x) const override
    {
        double scale = get_arg(0);
        return log(x) - 2 * log(scale) - 0.5 * pow(x / scale, 2);
    }

    double _logratio_(double xold, double xnew) const override
    {
        double scale = get_arg(0);
        return log(xnew / xold) + 0.5 * pow(scale, -2) * (pow(xold, 2) - pow(xnew, 2));
    }

    double _rand_(rnd_gen &prng) override
    {
        double scale = get_arg(0);
        return scale * sqrt(-2 * log(unif_gen(prng)));
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
    Exponential_Distribution() : distr_gen(1.0)
    {
        _init_defaults_("exp",
                        {"rate"}, {1.0},
                        0.0, NULL_DOUBLE, true);
    }
    ~Exponential_Distribution() {};

    void assign_reshaping(bool discrete, double lower, double upper) override
    {
        if ((!tools::is_na(lower)) && (lower <= 0))
            throw std::invalid_argument(
                "Invalid min value " + std::to_string(lower) + " for distribution \"" + get_name() + "\".");
        _assign_reshaping_(((discrete) ? "left" : NULL_STRING), lower, upper);
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

    virtual double _logratio_(double xold, double xnew) const override
    {
        double rate = get_arg(0);
        return rate * (xold - xnew);
    }

    virtual void _update_args_() override
    {
        distr_gen.param(std::exponential_distribution<double>::param_type(get_arg(0)));
    }

    virtual double _rand_(rnd_gen &prng) override
    {
        return distr_gen(prng);
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

const str_vec AVAILABLE_DISTRIBUTIONS = {
    "unif", "norm", "lnorm", "gamma", "weibull", "rayleigh", "exp"};

/* Auxiliary hidden template to return the chosen distribution. */
std::shared_ptr<Base_Distribution> distribution_switchboard(const std::string &name)
{
    const str_vec &options = AVAILABLE_DISTRIBUTIONS;

    if (!tools::is_na(name))
    {
        if (name == "unif")
            return std::make_shared<Uniform_Distribution>();

        else if (name == "norm")
            return std::make_shared<Normal_Distribution>();

        else if (name == "lnorm")
            return std::make_shared<Lognorm_Distribution>();

        else if (name == "gamma")
            return std::make_shared<Gamma_Distribution>();

        else if (name == "weibull")
            return std::make_shared<Weibull_Distribution>();

        else if (name == "rayleigh")
            return std::make_shared<Rayleigh_Distribution>();

        else if (name == "exp")
            return std::make_shared<Exponential_Distribution>();

        else
            throw std::invalid_argument(
                "Unknown distribution \"" + name + "\"." + tools::find_closest(name, options));
    }

    return std::make_shared<Base_Distribution>();
}
