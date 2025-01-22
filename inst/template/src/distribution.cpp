/*
 * distribution.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include <cmath>

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
        double sum = 1. / shape;
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
        double b = scaled_x + 1. - shape;
        double c = 1. / GAMMAP_FPMIN;
        double d = 1. / b;
        double h = d;
        double delta = 0.;

        for (int i = 1; ((fabs(delta - 1.) > GAMMAP_EPS) && (i < GAMMAP_ITER)); ++i)
        {
            double an = -i * (i - shape);
            b += 2.;
            d = an * d + b;
            if (fabs(d) < GAMMAP_FPMIN)
                d = GAMMAP_FPMIN;
            c = b + an / c;
            if (fabs(c) < GAMMAP_FPMIN)
                c = GAMMAP_FPMIN;
            d = 1. / d;
            delta = c * d;
            h *= delta;
        }
        return 1. - exp(-scaled_x + shape * log(scaled_x) - lgamma(shape)) * h;
    }

    /* Function to compute the gamma cumulative distribution function. */
    double gamma_cdf(double x, double shape, double scale)
    {
        double scaled_x = x / scale;
        if (scaled_x < shape + 10.)
            return gamma_cdf_series(scaled_x, shape);
        else
            return gamma_cdf_cf(scaled_x, shape);
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

double Distribution::pdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_pdf_)(x, opts);
}
double Distribution::cdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_cdf_)(x, opts);
}
double Distribution::logpdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_logpdf_)(x, opts);
}
double Distribution::logratio(double xold, double xnew, const dbl_vec &opts) const
{
    return (this->*_func_logratio_)(xold, xnew, opts);
}
double Distribution::rand(rnd_gen &prng)
{
    _update_args_();
    return (this->*_func_rand_)(prng);
}

/* Default reshaping initialization. */
void Distribution::init_reshaping(bool discrete)
{
    init_close_range();
    _init_discretization_((discrete) ? "center" : NULL_STRING);
    _init_truncation_((discrete) ? "center" : NULL_STRING);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution - virtual  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_pdf_(double x, const dbl_vec &opts) const
{
    throw std::invalid_argument(
        "Attempt to call undefined method \"pdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::_cdf_(double x, const dbl_vec &opts) const
{
    throw std::invalid_argument(
        "Attempt to call undefined method \"cdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::_logpdf_(double x, const dbl_vec &opts) const
{
    throw std::invalid_argument(
        "Attempt to call undefined method \"logpdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    UNUSED(xold);
    UNUSED(xnew);
    return 0.;
}
void Distribution::_update_args_() {};
double Distribution::_rand_(rnd_gen &prng)
{
    UNUSED(prng);
    throw std::invalid_argument(
        "Attempt to call undefined method \"rand\" of distribution \"" + show_name() + "\".");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : deaults  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_default_pdf_(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc_)(_pdf_(x, opts), opts);
}
double Distribution::_default_cdf_(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(x, opts), opts), opts);
}
double Distribution::_default_logpdf_(double x, const dbl_vec &opts) const
{
    return (this->*_func_logtrunc_)(_logpdf_(x, opts), opts);
}
double Distribution::_default_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    return _logratio_(xold, xnew, opts);
}

/* Check accept/reject within range is not exceeding maximum number of attempts. */
double Distribution::_default_rand_(rnd_gen &prng)
{
    double xnew = _rand_(prng);

    unsigned n_attempts = 0;
    while (!valid(xnew))
    {
        if (n_attempts > MAX_NUM_ATTEMPTS)
            throw std::runtime_error(
                "Exceeded attempts (N > " + std::to_string(MAX_NUM_ATTEMPTS) + ")" + " in random draw of " + show_name() + "(" + show_arg_values() + ")" + " with truncation {" + show_min() + ", " + show_max() + "}.");
        xnew = _rand_(prng);
        ++n_attempts;
    }

    return xnew;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized center /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_discr_center_pdf_(double x, const dbl_vec &opts) const
{
    double d = round(x);
    return (this->*_func_trunc_)(_cdf_(d + .5, opts) - _cdf_(d - .5, opts), opts);
}

double Distribution::_discr_center_cdf_(double x, const dbl_vec &opts) const
{
    double d = round(x);
    return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(d + .5, opts), opts), opts);
}

double Distribution::_discr_center_logpdf_(double x, const dbl_vec &opts) const
{
    return log(_discr_center_pdf_(x, opts));
}

double Distribution::_discr_center_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    double dold = round(xold), dnew = round(xnew);
    return log((_cdf_(dnew + .5, opts) - _cdf_(dnew - .5, opts)) / (_cdf_(dold + .5, opts) - _cdf_(dold - .5, opts)));
}

double Distribution::_discr_center_rand_(rnd_gen &prng)
{
    // UPDATE : lround to avoid casting like round(-0.1) -> -0.0
    return static_cast<double>(std::lround(_default_rand_(prng)));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : discretized left /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_discr_left_pdf_(double x, const dbl_vec &opts) const
{
    double d = floor(x);
    return (this->*_func_trunc_)(_cdf_(d + 1., opts) - _cdf_(d, opts), opts);
}

double Distribution::_discr_left_cdf_(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc_)((this->*_func_shift_)(_cdf_(floor(x) + 1., opts), opts), opts);
}

double Distribution::_discr_left_logpdf_(double x, const dbl_vec &opts) const
{
    return log(_discr_left_pdf_(x, opts));
}

double Distribution::_discr_left_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    double dold = floor(xold), dnew = floor(xnew);
    return log((_cdf_(dnew + 1., opts) - _cdf_(dnew, opts)) / (_cdf_(dold + 1., opts) - _cdf_(dold, opts)));
}

double Distribution::_discr_left_rand_(rnd_gen &prng) { return floor(_default_rand_(prng)); }

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Base Distribution : truncation (functional)  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_trunc_dummy_(double x, const dbl_vec &opts) const
{
    UNUSED(opts);
    return x;
}
double Distribution::_trunc_min_(double x, const dbl_vec &opts) const
{
    return x / (1. - _cdf_(get_min(), opts));
}
double Distribution::_trunc_max_(double x, const dbl_vec &opts) const
{
    return x / (_cdf_(get_max(), opts));
}
double Distribution::_trunc_both_sides_(double x, const dbl_vec &opts) const
{
    return x / (_cdf_(get_max(), opts) - _cdf_(get_max(), opts));
}

// case of center discretization
double Distribution::_trunc_min_discr_center_(double x, const dbl_vec &opts) const
{
    return x / (1. - _cdf_(floor(get_min()) - .5, opts));
}
double Distribution::_trunc_max_discr_center_(double x, const dbl_vec &opts) const
{
    return x / (_cdf_(ceil(get_max()) + .5, opts));
}
double Distribution::_trunc_both_sides_discr_center_(double x, const dbl_vec &opts) const
{
    return x / (_cdf_(ceil(get_max()) + .5, opts) - _cdf_(floor(get_min()) - .5, opts));
}

// case of left discretization
double Distribution::_trunc_min_discr_left_(double x, const dbl_vec &opts) const
{
    return x / (1. - _cdf_(floor(get_min()), opts));
}
double Distribution::_trunc_max_discr_left_(double x, const dbl_vec &opts) const
{
    return x / (_cdf_(ceil(get_max()), opts));
}
double Distribution::_trunc_both_sides_discr_left_(double x, const dbl_vec &opts) const
{
    return x / (_cdf_(ceil(get_max()), opts) - _cdf_(floor(get_min()), opts));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Base Distribution : logtruncation (functional)  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_logtrunc_dummy_(double x, const dbl_vec &opts) const
{
    UNUSED(opts);
    return x;
}
double Distribution::_logtrunc_min_(double x, const dbl_vec &opts) const { return x - log(1. - _cdf_(get_min(), opts)); }
double Distribution::_logtrunc_max_(double x, const dbl_vec &opts) const { return x - log(_cdf_(get_max(), opts)); }
double Distribution::_logtrunc_both_sides_(double x, const dbl_vec &opts) const { return x - log(_cdf_(get_max(), opts) - _cdf_(get_min(), opts)); }

// case of center discretization
double Distribution::_logtrunc_min_discr_center_(double x, const dbl_vec &opts) const { return x - log(1. - _cdf_(floor(get_min()) - .5, opts)); }
double Distribution::_logtrunc_max_discr_center_(double x, const dbl_vec &opts) const { return x - log(_cdf_(ceil(get_max()) + .5, opts)); }
double Distribution::_logtrunc_both_sides_discr_center_(double x, const dbl_vec &opts) const { return x - log(_cdf_(ceil(get_max()) + .5, opts) - _cdf_(floor(get_min()) - .5, opts)); }

// case of left discretization
double Distribution::_logtrunc_min_discr_left_(double x, const dbl_vec &opts) const { return x - log(1. - _cdf_(floor(get_min()), opts)); }
double Distribution::_logtrunc_max_discr_left_(double x, const dbl_vec &opts) const { return x - log(_cdf_(ceil(get_max()), opts)); }
double Distribution::_logtrunc_both_sides_discr_left_(double x, const dbl_vec &opts) const { return x - log(_cdf_(ceil(get_max()), opts) - _cdf_(floor(get_min()), opts)); }

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  Base Distribution : shift cumulative (functional)  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_shift_dummy_(double x, const dbl_vec &opts) const { return x; }
double Distribution::_shift_min_(double x, const dbl_vec &opts) const { return x - _cdf_(get_min(), opts); }
double Distribution::_shift_min_discr_center_(double x, const dbl_vec &opts) const { return x - _cdf_(floor(get_min()) - .5, opts); }
double Distribution::_shift_min_discr_left_(double x, const dbl_vec &opts) const { return x - _cdf_(floor(get_min()), opts); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : assign reshaping  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Select functional method in case of discretization. */
void Distribution::_init_discretization_(const std::string &discrete)
{
    //===== DISCRETIZATION =====//
    if (tools::is_na(discrete)) // no discretization (default)
    {
        _func_pdf_ = &Distribution::_default_pdf_;
        _func_cdf_ = &Distribution::_default_cdf_;
        _func_logpdf_ = &Distribution::_default_logpdf_;
        _func_logratio_ = &Distribution::_default_logratio_;
        _func_rand_ = &Distribution::_default_rand_;
    }
    else if (discrete == "center") // centered discretization
    {
        _func_pdf_ = &Distribution::_discr_center_pdf_;
        _func_cdf_ = &Distribution::_discr_center_cdf_;
        _func_logpdf_ = &Distribution::_discr_center_logpdf_;
        _func_logratio_ = &Distribution::_discr_center_logratio_;
        _func_rand_ = &Distribution::_discr_center_rand_;
    }
    else if (discrete == "left") // left adjusted discretization
    {
        _func_pdf_ = &Distribution::_discr_left_pdf_;
        _func_cdf_ = &Distribution::_discr_left_cdf_;
        _func_logpdf_ = &Distribution::_discr_left_logpdf_;
        _func_logratio_ = &Distribution::_discr_left_logratio_;
        _func_rand_ = &Distribution::_discr_left_rand_;
    }
    else
    {
        throw std::invalid_argument(
            "Unrecognized discretization method \"" + discrete + "\".");
    }
}

/* Select functional method in case of truncation. */
void Distribution::_init_truncation_(const std::string &discrete)
{
    //----- TRUNCATION -----//
    switch (which_valid_case())
    {
    case 1:                       // lower range truncation
        if (discrete == "center") // "center" discretization
        {
            _func_trunc_ = &Distribution::_trunc_min_discr_center_;
            _func_logtrunc_ = &Distribution::_logtrunc_min_discr_center_;
            _func_shift_ = &Distribution::_shift_min_discr_center_;
        }
        else if (discrete == "left") // "left" discretization
        {
            _func_trunc_ = &Distribution::_trunc_min_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_min_discr_left_;
            _func_shift_ = &Distribution::_shift_min_discr_left_;
        }
        else // no discretization
        {
            _func_trunc_ = &Distribution::_trunc_min_;
            _func_logtrunc_ = &Distribution::_logtrunc_min_;
            _func_shift_ = &Distribution::_shift_min_;
        }
        break;
    case 2:                       // upper range truncation
        if (discrete == "center") // "center" discretization
        {
            _func_trunc_ = &Distribution::_trunc_max_discr_center_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_discr_center_;
        }
        else if (discrete == "left") // "left" discretization
        {
            _func_trunc_ = &Distribution::_trunc_max_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_discr_left_;
        }
        else // no discretization
        {
            _func_trunc_ = &Distribution::_trunc_max_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_;
        }
        _func_shift_ = &Distribution::_shift_dummy_;
        break;
    case 3:                       // both sides truncation
        if (discrete == "center") // "center" discretization
        {
            _func_trunc_ = &Distribution::_trunc_both_sides_discr_center_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_discr_center_;
            _func_shift_ = &Distribution::_shift_min_discr_center_;
        }
        else if (discrete == "left") // "left" discretization
        {
            _func_trunc_ = &Distribution::_trunc_both_sides_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_discr_left_;
            _func_shift_ = &Distribution::_shift_min_discr_left_;
        }
        else // no discretization
        {
            _func_trunc_ = &Distribution::_trunc_both_sides_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_;
            _func_shift_ = &Distribution::_shift_min_;
        }
        break;
    default: // case 0: no truncation
        _func_trunc_ = &Distribution::_trunc_dummy_;
        _func_logtrunc_ = &Distribution::_logtrunc_dummy_;
        _func_shift_ = &Distribution::_shift_dummy_;
        break;
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

class Uniform_Distribution : public Distribution
{
private:
    std::uniform_real_distribution<double> distr_gen;

public:
    //===== CONSTRUCTOR =====//
    Uniform_Distribution() : distr_gen(0., 1.)
    {
        _init_defaults_("unif", {"min", "max"}, {0., 1.},
                        NULL_DOUBLE, NULL_DOUBLE, false);
    }
    ~Uniform_Distribution() = default;

protected:
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(x);
        UNUSED(opts);
        return 1.;
    }

    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return x;
    }

    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(x);
        UNUSED(opts);
        return 0.;
    }

    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(xold);
        UNUSED(xnew);
        UNUSED(opts);
        return 0.;
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
class Normal_Distribution : public Distribution
{
private:
    std::normal_distribution<double> distr_gen;

public:
    //===== CONSTRUCTOR =====//
    Normal_Distribution() : distr_gen(0., 1.)
    {
        _init_defaults_("norm", {"mean", "sd"}, {0., 1.},
                        NULL_DOUBLE, NULL_DOUBLE, false);
    }
    ~Normal_Distribution() = default;

protected:
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        UNUSED(opts);
        return .39894228 * exp(-.5 * pow((x - mean) / sd, 2)) / sd;
    }

    double _cdf_(double x, const dbl_vec &opts) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        UNUSED(opts);
        return .5 + .5 * erf(.707106781 * (x - mean) / sd);
    }

    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        UNUSED(opts);
        return -.5 * pow((x - mean) / sd, 2) - .918938533 - log(sd);
    }

    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        double mean = get_arg(0);
        double sd = get_arg(1);
        UNUSED(opts);
        return .5 * (pow((xold - mean) / sd, 2) - pow((xnew - mean) / sd, 2));
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

class Lognorm_Distribution : public Distribution
{
private:
    std::lognormal_distribution<double> distr_gen;

public:
    //===== CONSTRUCTOR =====//
    Lognorm_Distribution() : distr_gen(0., 1.)
    {
        _init_defaults_("lnorm", {"meanlog", "sdlog"}, {0., 1.},
                        0., NULL_DOUBLE, true);
    }
    ~Lognorm_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? "left" : NULL_STRING));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() == 0) && (tools::is_na(get_max())))
            init_open_range();
        else
        {
            init_close_range();
            _init_truncation_(((discrete) ? "left" : NULL_STRING));
        }
    }

protected:
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        UNUSED(opts);
        return exp(-.5 * pow((log(x) - meanlog) / sdlog, 2)) * .39894228 / (x * sdlog);
    }

    double _cdf_(double x, const dbl_vec &opts) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        UNUSED(opts);
        return .5 + .5 * erf(.707106781 * (log(x) - meanlog) / sdlog);
    }

    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        UNUSED(opts);
        return -.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x * sdlog) - .918938533;
    }

    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        double meanlog = get_arg(0);
        double sdlog = get_arg(1);
        UNUSED(opts);
        return .5 * (pow((log(xold) - meanlog) / sdlog, 2) - pow((log(xnew) - meanlog) / sdlog, 2)) + log(xold / xnew);
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

class Gamma_Distribution : public Distribution
{
private:
    std::gamma_distribution<double> distr_gen;

public:
    //===== CONSTRUCTOR =====//
    Gamma_Distribution() : distr_gen(1., 1.)
    {
        _init_defaults_("gamma", {"shape", "scale"}, {1., 1.},
                        0., NULL_DOUBLE, true);
    }
    ~Gamma_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? "left" : NULL_STRING));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() == 0) && (tools::is_na(get_max())))
            init_close_range();
        else
        {
            init_close_range();
            _init_truncation_(((discrete) ? "left" : NULL_STRING));
        }
    }

protected:
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        UNUSED(opts);
        return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
    }

    double _cdf_(double x, const dbl_vec &opts) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        UNUSED(opts);
        return auxgammap::gamma_cdf(x, shape, scale);
    }

    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        UNUSED(opts);
        return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
    }

    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        UNUSED(opts);
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

class Weibull_Distribution : public Distribution
{
private:
    std::weibull_distribution<double> distr_gen;

public:
    //===== CONSTRUCTOR =====//
    Weibull_Distribution() : distr_gen(1., 1.)
    {
        _init_defaults_("weibull", {"shape", "scale"}, {1., 1.},
                        0., NULL_DOUBLE, true);
    }
    ~Weibull_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? "left" : NULL_STRING));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else
        {
            init_close_range();
            if ((get_min() > 0) || (!tools::is_na(get_max())))
                _init_truncation_(((discrete) ? "left" : NULL_STRING));
        }
    }

protected:
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        UNUSED(opts);
        return (shape / scale) * pow((x / scale), shape - 1) * exp(-pow(x / scale, shape));
    }

    double _cdf_(double x, const dbl_vec &opts) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        UNUSED(opts);
        return 1. - exp(-pow(x / scale, shape));
    }

    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        UNUSED(opts);
        return log(shape / scale) + (shape - 1) * log(x / scale) - pow(x / scale, shape);
    }

    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        double shape = get_arg(0);
        double scale = get_arg(1);
        UNUSED(opts);
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

class Rayleigh_Distribution : public Distribution
{
private:
    std::uniform_real_distribution<double> unif_gen;

public:
    //===== CONSTRUCTOR =====//
    Rayleigh_Distribution() : unif_gen(0., 1.)
    {
        _init_defaults_("rayleigh", {"scale"}, {1.},
                        0., NULL_DOUBLE, true);
    }
    ~Rayleigh_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? "left" : NULL_STRING));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else
        {
            init_close_range();
            if ((get_min() > 0) || (!tools::is_na(get_max())))
                _init_truncation_(((discrete) ? "left" : NULL_STRING));
        }
    }

protected:
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        double scale = get_arg(0);
        UNUSED(opts);
        return x * pow(scale, -2) * exp(-.5 * pow(x / scale, 2));
    }

    double _cdf_(double x, const dbl_vec &opts) const override
    {
        double scale = get_arg(0);
        UNUSED(opts);
        return 1. - exp(-.5 * pow(x / scale, 2));
    }

    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        double scale = get_arg(0);
        UNUSED(opts);
        return log(x) - 2 * log(scale) - .5 * pow(x / scale, 2);
    }

    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        double scale = get_arg(0);
        UNUSED(opts);
        return log(xnew / xold) + .5 * pow(scale, -2) * (pow(xold, 2) - pow(xnew, 2));
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

class Exponential_Distribution : public Distribution
{
private:
    std::exponential_distribution<double> distr_gen;

public:
    //===== CONSTRUCTOR =====//
    Exponential_Distribution() : distr_gen(1.)
    {
        _init_defaults_("exp", {"rate"}, {1.},
                        0., NULL_DOUBLE, true);
    }
    ~Exponential_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? "left" : NULL_STRING));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else
        {
            init_close_range();
            if ((get_min() > 0) || (!tools::is_na(get_max())))
                _init_truncation_(((discrete) ? "left" : NULL_STRING));
        }
    }

protected:
    virtual double _pdf_(double x, const dbl_vec &opts) const override
    {
        double rate = get_arg(0);
        UNUSED(opts);
        return rate * exp(-rate * x);
    }

    virtual double _cdf_(double x, const dbl_vec &opts) const override
    {
        double rate = get_arg(0);
        UNUSED(opts);
        return 1. - exp(-rate * x);
    }

    virtual double _logpdf_(double x, const dbl_vec &opts) const override
    {
        double rate = get_arg(0);
        UNUSED(opts);
        return log(rate) - rate * x;
    }

    virtual double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        double rate = get_arg(0);
        UNUSED(opts);
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

/* Auxiliary hidden template to return the chosen distribution. */
std::shared_ptr<Distribution> distribution_switchboard(const std::string &name)
{
    const str_vec &options = {
        "unif", "norm", "lnorm", "gamma", "weibull", "rayleigh", "exp"};

    if (tools::is_na(name))
        return std::make_shared<Distribution>();

    else if (name == "unif")
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
