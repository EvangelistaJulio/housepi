/*
 * distribution.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include <cmath>

#include "binomial.h"
#include "distribution.h"

#define D_DISCR_TYPE_NONE 0 // none
#define D_DISCR_TYPE_LEFT 1 // left
#define D_DISCR_TYPE_CNTR 2 // center
#define D_DISCR_TYPE_DFND 3 // defined

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
    init_range(false, false);
    _init_discretization_((discrete) ? D_DISCR_TYPE_CNTR : D_DISCR_TYPE_NONE);
    _init_truncation_((discrete) ? D_DISCR_TYPE_CNTR : D_DISCR_TYPE_NONE);
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
///////////////////////////  Base Distribution : discretized defined ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_discr_defined_pdf_(double x, const dbl_vec &opts) const
{
    return _default_pdf_(floor(x), opts);
}
double Distribution::_discr_defined_cdf_(double x, const dbl_vec &opts) const
{
    return _default_pdf_(floor(x), opts);
}
double Distribution::_discr_defined_logpdf_(double x, const dbl_vec &opts) const
{
    return _default_logpdf_(floor(x), opts);
}
double Distribution::_discr_defined_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    return _default_logratio_(floor(xold), floor(xnew), opts);
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

/* Tell if the distribution is discrete. */
bool Distribution::is_discrete() const
{
    return _func_pdf_ == &Distribution::_default_pdf_;
}

/* Select functional method in case of discretization. */
void Distribution::_init_discretization_(ushort discr_type)
{
    //===== DISCRETIZATION =====//
    switch (discr_type)
    {
    case D_DISCR_TYPE_NONE: // no discretization (default)
        _func_pdf_ = &Distribution::_default_pdf_;
        _func_cdf_ = &Distribution::_default_cdf_;
        _func_logpdf_ = &Distribution::_default_logpdf_;
        _func_logratio_ = &Distribution::_default_logratio_;
        _func_rand_ = &Distribution::_default_rand_;
        break;

    case D_DISCR_TYPE_DFND: // discrete defined
        _func_pdf_ = &Distribution::_discr_defined_pdf_;
        _func_cdf_ = &Distribution::_discr_defined_cdf_;
        _func_logpdf_ = &Distribution::_discr_defined_logpdf_;
        _func_logratio_ = &Distribution::_discr_defined_logratio_;
        _func_rand_ = &Distribution::_default_rand_;
        break;

    case D_DISCR_TYPE_CNTR: // centered discretization
        _func_pdf_ = &Distribution::_discr_center_pdf_;
        _func_cdf_ = &Distribution::_discr_center_cdf_;
        _func_logpdf_ = &Distribution::_discr_center_logpdf_;
        _func_logratio_ = &Distribution::_discr_center_logratio_;
        _func_rand_ = &Distribution::_discr_center_rand_;
        break;

    case D_DISCR_TYPE_LEFT: // left discretization
        _func_pdf_ = &Distribution::_discr_left_pdf_;
        _func_cdf_ = &Distribution::_discr_left_cdf_;
        _func_logpdf_ = &Distribution::_discr_left_logpdf_;
        _func_logratio_ = &Distribution::_discr_left_logratio_;
        _func_rand_ = &Distribution::_discr_left_rand_;
        break;
    }
}

/* Select functional method in case of truncation. */
void Distribution::_init_truncation_(ushort discr_type)
{
    //----- TRUNCATION -----//
    switch (which_valid_case())
    {
    case 1: //----- lower range truncation -----//
        switch (discr_type)
        {
        case D_DISCR_TYPE_CNTR: // centered discretization
            _func_trunc_ = &Distribution::_trunc_min_discr_center_;
            _func_logtrunc_ = &Distribution::_logtrunc_min_discr_center_;
            _func_shift_ = &Distribution::_shift_min_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc_ = &Distribution::_trunc_min_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_min_discr_left_;
            _func_shift_ = &Distribution::_shift_min_discr_left_;
            break;
        default: // no discretization or D_DISCR_TYPE_DFND
            _func_trunc_ = &Distribution::_trunc_min_;
            _func_logtrunc_ = &Distribution::_logtrunc_min_;
            _func_shift_ = &Distribution::_shift_min_;
        }
        break;
    case 2: //----- upper range truncation -----//
        switch (discr_type)
        {
        case D_DISCR_TYPE_CNTR: // centered discretization
            _func_trunc_ = &Distribution::_trunc_max_discr_center_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc_ = &Distribution::_trunc_max_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_discr_left_;
            break;
        default: // no discretization or D_DISCR_TYPE_DFND
            _func_trunc_ = &Distribution::_trunc_max_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_;
        }
        _func_shift_ = &Distribution::_shift_dummy_;
        break;
    case 3: //----- both sides truncation -----//
        switch (discr_type)
        {
        case D_DISCR_TYPE_CNTR: // centered discretization
            _func_trunc_ = &Distribution::_trunc_both_sides_discr_center_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_discr_center_;
            _func_shift_ = &Distribution::_shift_min_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc_ = &Distribution::_trunc_both_sides_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_discr_left_;
            _func_shift_ = &Distribution::_shift_min_discr_left_;
            break;
        default: // no discretization or D_DISCR_TYPE_DFND
            _func_trunc_ = &Distribution::_trunc_both_sides_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_;
            _func_shift_ = &Distribution::_shift_min_;
        }
        break;
    default: //----- no truncation (case 0) -----//
        _func_trunc_ = &Distribution::_trunc_dummy_;
        _func_logtrunc_ = &Distribution::_logtrunc_dummy_;
        _func_shift_ = &Distribution::_shift_dummy_;
        break;
    }
}

//===== SHOW =====//
const std::map<double, double> Distribution::showcase_function(const dbl_vec &opts) const
{
    std::map<double, double> output = {};

    // pick extremes for the showcase
    bool is_min_open = false;
    bool is_max_open = false;
    double lower = this->get_min();
    double upper = this->get_max();

    // adjust extremes if unbound
    if ((tools::is_na(lower)) && (tools::is_na(upper)))
    {
        upper = 10.;
        lower = -upper;
    }
    else
    {
        // check if extrems are open
        is_min_open = !this->valid(lower);
        is_max_open = !this->valid(upper);
    }

    double delta_x = (upper - lower) * 1.e-10;

    // adjust extremes if interval is open
    if (is_min_open)
        lower += delta_x;
    if (is_max_open)
        upper -= delta_x;

    // compute values to showcase
    for (double x = lower; x <= upper; x += delta_x)
        output.emplace(x, this->pdf(x, opts));

    return output;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.3989422804014327 /* 1/sqrt(2*pi) */
#endif

#ifndef M_1_SQRT_2
#define M_1_SQRT_2 0.7071067811865476 /* 1/sqrt(2) */
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.9189385332046727 /* log(sqrt(2*pi)) */
#endif

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Uniform Distribution : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Uniform Distribution derived class. */
class Uniform_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Uniform_Distribution() : distr_gen(0., 1.)
    {
        _init_defaults_("unif", {"min", "max"}, {0., 1.},
                        NULL_DOUBLE, NULL_DOUBLE, false, false);
    }
    ~Uniform_Distribution() = default;

protected:
    std::uniform_real_distribution<double> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::uniform_real_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
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
    double _rand_(rnd_gen &prng) override
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
///////////////////////////////  Normal Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace normal_impl
{
    inline double pdf(double x, double mean, double sd)
    {
        return M_1_SQRT_2PI * exp(-.5 * pow((x - mean) / sd, 2)) / sd;
    }
    inline double cdf(double x, double mean, double sd)
    {
        return .5 + .5 * erf(M_1_SQRT_2 * (x - mean) / sd);
    }
    inline double logpdf(double x, double mean, double sd)
    {
        return -.5 * pow((x - mean) / sd, 2) - M_LN_SQRT_2PI - log(sd);
    }
    inline double logratio(double xold, double xnew, double mean, double sd)
    {
        return .5 * (pow((xold - mean) / sd, 2) - pow((xnew - mean) / sd, 2));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Normal Distribution : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Normal Distribution derived class. */
class Normal_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Normal_Distribution() : distr_gen(0., 1.)
    {
        _init_defaults_("norm", {"mean", "sd"}, {0., 1.},
                        NULL_DOUBLE, NULL_DOUBLE, false, false);
    }
    ~Normal_Distribution() = default;

protected:
    std::normal_distribution<double> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::normal_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return normal_impl::pdf(x, get_arg(0), get_arg(1));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return normal_impl::cdf(x, get_arg(0), get_arg(1));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return normal_impl::logpdf(x, get_arg(0), get_arg(1));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return normal_impl::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double _rand_(rnd_gen &prng) override
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
//////////////////////////////  Lognormal Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace lognormal_impl
{
    inline double pdf(double x, double meanlog, double sdlog)
    {
        return exp(-.5 * pow((log(x) - meanlog) / sdlog, 2)) * M_1_SQRT_2PI / (x * sdlog);
    }
    inline double cdf(double x, double meanlog, double sdlog)
    {
        return .5 + .5 * erf(M_1_SQRT_2 * (log(x) - meanlog) / sdlog);
    }
    inline double logpdf(double x, double meanlog, double sdlog)
    {
        return -.5 * pow((log(x) - meanlog) / sdlog, 2) - log(x * sdlog) - M_LN_SQRT_2PI;
    }
    inline double logratio(double xold, double xnew, double meanlog, double sdlog)
    {
        return .5 * (pow((log(xold) - meanlog) / sdlog, 2) - pow((log(xnew) - meanlog) / sdlog, 2)) + log(xold / xnew);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Lognormal Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Lognormal Distribution derived class. */
class Lognormal_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Lognormal_Distribution() : distr_gen(0., 1.)
    {
        _init_defaults_("lnorm", {"meanlog", "sdlog"}, {0., 1.},
                        0., NULL_DOUBLE, true, false);
    }
    ~Lognormal_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0) || (!tools::is_na(get_max())))
        {
            init_range((get_min() == 0), false);
            _init_truncation_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));
        }
    }

protected:
    std::lognormal_distribution<double> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::lognormal_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return lognormal_impl::pdf(x, get_arg(0), get_arg(1));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return lognormal_impl::cdf(x, get_arg(0), get_arg(1));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return lognormal_impl::logpdf(x, get_arg(0), get_arg(1));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return lognormal_impl::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double _rand_(rnd_gen &prng) override
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
////////////////////////  pgamma : regularized incomplete gamma function  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace gamma_impl_core
{
    /*
     * FIXME: To better reference this:
     * Press, Teukolsky, Vetterling, and Flannery
     * Numerical Recipes 3rd Edition: The Art of Scientific Computing
     */
    constexpr int MAX_ITER = 1000;
    constexpr double EPSILON = 1.e-12;
    constexpr double MIN_DBL = 1.e-18;

    /* Lower regulairzed incomplete gamma function using series expansion. */
    inline double _pgamma_series_(double scaled_x, double shape)
    {
        double sum = 1. / shape;
        double term = sum;
        for (int i = 1; (term > EPSILON * sum) && (i <= MAX_ITER); ++i)
        {
            term *= scaled_x / (shape + i);
            sum += term;
        }
        return sum * exp(-scaled_x + shape * log(scaled_x) - lgamma(shape));
    }

    /* Upper regularized incomplete gamma function using continued fraction. */
    inline double _pgamma_continued_fraction_(double scaled_x, double shape)
    {
        double b = scaled_x + 1. - shape;
        double c = 1. / MIN_DBL;
        double d = 1. / b;
        double h = d;
        double delta = 0.;

        for (int i = 1; ((fabs(delta - 1.) > EPSILON) && (i <= MAX_ITER)); ++i)
        {
            double an = -i * (i - shape);
            b += 2.;
            d = an * d + b;
            if (fabs(d) < MIN_DBL)
                d = MIN_DBL;
            c = b + an / c;
            if (fabs(c) < MIN_DBL)
                c = MIN_DBL;
            d = 1. / d;
            delta = c * d;
            h *= delta;
        }

        return 1. - exp(-scaled_x + shape * log(scaled_x) - lgamma(shape)) * h;
    }

    /* Regularized incomplete gamma function. */
    inline double pgamma(double x, double shape)
    {
        if (x < shape + 10.)
            return _pgamma_series_(x, shape);
        else
            return _pgamma_continued_fraction_(x, shape);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Gamma Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace gamma_impl
{
    inline double pdf(double x, double shape, double scale)
    {
        return exp(-x / scale) * pow(x, shape - 1.) / (tgamma(shape) * pow(scale, shape));
    }
    inline double logpdf(double x, double shape, double scale)
    {
        return (shape - 1) * log(x) - lgamma(shape) - shape * log(scale) - x / scale;
    }
    inline double logratio(double xold, double xnew, double shape, double scale)
    {
        return (shape - 1) * (log(xnew / xold)) + (xold - xnew) / scale;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Gamma Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Gamma Distribution derived class. */
class Gamma_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Gamma_Distribution() : distr_gen(1., 1.)
    {
        _init_defaults_("gamma", {"shape", "scale"}, {1., 1.},
                        0., NULL_DOUBLE, false, false);
    }
    ~Gamma_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0) || (!tools::is_na(get_max())))
        {
            init_range(false, false);
            _init_truncation_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));
        }
    }

protected:
    std::gamma_distribution<double> distr_gen;

    //==== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::gamma_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //==== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return gamma_impl::pdf(x, get_arg(0), get_arg(1));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return gamma_impl_core::pgamma(x / get_arg(1), get_arg(0));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return gamma_impl::logpdf(x, get_arg(0), get_arg(1));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return gamma_impl::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double _rand_(rnd_gen &prng) override
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
///////////////////////////////  Weibull Distribution : methods  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace weibull_impl
{
    inline double _pdf_(double x, double shape, double scale)
    {
        return (shape / scale) * pow((x / scale), shape - 1.) * exp(-pow(x / scale, shape));
    }
    inline double _cdf_(double x, double shape, double scale)
    {
        return 1. - exp(-pow(x / scale, shape));
    }
    inline double _logpdf_(double x, double shape, double scale)
    {
        return log(shape / scale) + (shape - 1.) * log(x / scale) - pow(x / scale, shape);
    }
    inline double _logratio_(double xold, double xnew, double shape, double scale)
    {
        return (shape - 1.) * (log(xnew / xold)) + pow(xold / scale, shape) - pow(xnew / scale, shape);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Weibull Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Weibull Distribution derived class. */
class Weibull_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Weibull_Distribution() : distr_gen(1., 1.)
    {
        _init_defaults_("weibull", {"shape", "scale"}, {1., 1.},
                        0., NULL_DOUBLE, false, false);
    }
    ~Weibull_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0) || (!tools::is_na(get_max())))
        {
            init_range(false, false);
            _init_truncation_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));
        }
    }

protected:
    std::weibull_distribution<double> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::weibull_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return weibull_impl::_pdf_(x, get_arg(0), get_arg(1));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return weibull_impl::_cdf_(x, get_arg(0), get_arg(1));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return weibull_impl::_logpdf_(x, get_arg(0), get_arg(1));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return weibull_impl::_logratio_(xold, xnew, get_arg(0), get_arg(1));
    }
    double _rand_(rnd_gen &prng) override
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
///////////////////////////////  Rayleigh Distribution : methods  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace rayleigh_impl
{
    inline double pdf(double x, double scale)
    {
        return x * pow(scale, -2) * exp(-.5 * pow(x / scale, 2));
    }
    inline double cdf(double x, double scale)
    {
        return 1. - exp(-.5 * pow(x / scale, 2));
    }
    inline double logpdf(double x, double scale)
    {
        return log(x) - 2 * log(scale) - .5 * pow(x / scale, 2);
    }
    inline double logratio(double xold, double xnew, double scale)
    {
        return log(xnew / xold) + .5 * pow(scale, -2) * (pow(xold, 2) - pow(xnew, 2));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Rayleigh Distribution : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Rayleigh Distribution derived class. */
class Rayleigh_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Rayleigh_Distribution() : unif_gen(0., 1.)
    {
        _init_defaults_("rayleigh", {"scale"}, {1.},
                        0., NULL_DOUBLE, false, false);
    }
    ~Rayleigh_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0) || (!tools::is_na(get_max())))
        {
            init_range(false, false);
            _init_truncation_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));
        }
    }

protected:
    std::uniform_real_distribution<double> unif_gen;

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return rayleigh_impl::pdf(x, get_arg(0));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return rayleigh_impl::cdf(x, get_arg(0));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return rayleigh_impl::logpdf(x, get_arg(0));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return rayleigh_impl::logratio(xold, xnew, get_arg(0));
    }
    double _rand_(rnd_gen &prng) override
    {
        return get_arg(0) * sqrt(-2 * log(unif_gen(prng)));
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  ibeta : regularized incomplete beta function  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace beta_impl_core
{
    inline constexpr double EPSILON = 1e-12;
    inline constexpr int MAX_ITER = 1000;

    /* logbeta function. */
    inline double _lbeta_(double a, double b)
    {
        return lgamma(a) + lgamma(b) - lgamma(a + b);
    }

    /* Beta function. */
    inline double _beta_(double a, double b)
    {
        return exp(_lbeta_(a, b));
    }

    /* Continued fraction evaluation for regularized incomplete beta. */
    inline double _ibeta_continued_fraction_(double x, double a, double b)
    {
        double qab = a + b;
        double qap = a + 1.;
        double qam = a - 1.;
        double c = 1.;
        double d = 1. / (1. - qab * x / qap);
        double h = d;

        for (unsigned m = 1; m <= MAX_ITER; ++m)
        {
            double m2 = 2. * m;
            double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1. / (1. + aa * d);
            c = 1. + aa / c;
            h *= d * c;
            if (fabs(d - 1.) < EPSILON)
                return h;
        }

        return std::numeric_limits<double>::quiet_NaN();
    }

    /* Series expansion for small x. */
    inline double _ibeta_series_(double x, double a, double b)
    {
        double sum = 1.;
        double term = 1.;
        double apb = a + b;
        double ap = a + 1.;

        for (unsigned n = 1; n <= MAX_ITER; ++n)
        {
            term *= (apb * x) / (ap * n);
            sum += term;
            if (fabs(term) < fabs(sum) * EPSILON)
                break;
            ++ap;
        }

        return pow(x, a) * pow(1. - x, b) / (a * _beta_(a, b)) * sum;
    }

    // Regularized incomplete beta function (normalized)
    inline double ibeta(double x, double a, double b)
    {
        /* FIXME: should I implement parameter validation somewhere ?
        if (a <= 0 || b <= 0 || x < 0 || x > 1)
            return std::numeric_limits<double>::quiet_NaN();
        */

        // edge cases
        if (x == 0.)
            return 0.;
        if (x == 1.)
            return 1.;

        // choose computation method based on x
        if (x < (a + 1.) / (a + b + 2.))
            return _ibeta_series_(x, a, b);
        else
            return _ibeta_continued_fraction_(x, a, b);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Beta Distribution : methods  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace beta_impl
{
    inline double pdf(double x, double shape1, double shape2)
    {
        return exp((shape1 - 1.) * log(x) + (shape2 - 1.) * log(1. - x) - beta_impl_core::_lbeta_(shape1, shape2));
    }
    inline double logpdf(double x, double shape1, double shape2)
    {
        return (shape1 - 1.) * log(x) + (shape2 - 1.) * log(1. - x) - beta_impl_core::_lbeta_(shape1, shape2);
    }
    inline double logratio(double xold, double xnew, double shape1, double shape2)
    {
        return (shape1 - 1.) * log(xnew / xold) + (shape2 - 1.) * log((1. - xnew) / (1. - xold));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Beta Distribution : class  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Beta Distribution derived class. */
class Beta_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Beta_Distribution() : gamma_1(1., 1.), gamma_2(1., 1.)
    {
        _init_defaults_("beta", {"shape1", "shape2"}, {1., 1.},
                        0., 1., false, false);
    }
    ~Beta_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        if (discrete)
            throw std::invalid_argument(
                "Invalid option discrete for distribution \"" + show_name() + "\".");
        _init_discretization_(D_DISCR_TYPE_NONE);

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if (get_max() > 1)
            throw std::invalid_argument(
                "Invalid max value " + show_max() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0) || (get_max() < 1))
        {
            init_range(false, false);
            _init_truncation_(D_DISCR_TYPE_NONE);
        }
    }

protected:
    std::gamma_distribution<double> gamma_1;
    std::gamma_distribution<double> gamma_2;

    //==== SETTER =====//
    void _update_args_() override
    {
        gamma_1.param(std::gamma_distribution<double>::param_type(get_arg(0), 1.));
        gamma_2.param(std::gamma_distribution<double>::param_type(get_arg(1), 1.));
    }

    //==== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return beta_impl::pdf(x, get_arg(0), get_arg(1));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return beta_impl_core::ibeta(x, get_arg(0), get_arg(1));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return beta_impl::logpdf(x, get_arg(0), get_arg(1));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return beta_impl::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double _rand_(rnd_gen &prng) override
    {
        double x1 = gamma_1(prng);
        double x2 = gamma_2(prng);
        return x1 / (x1 + x2);
    }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Exponential  Distribution : methods  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace exponential_impl
{
    inline double pdf(double x, double rate)
    {
        return rate * exp(-rate * x);
    }
    inline double cdf(double x, double rate)
    {
        return 1. - exp(-rate * x);
    }
    inline double logpdf(double x, double rate)
    {
        return log(rate) - rate * x;
    }
    inline double logratio(double xold, double xnew, double rate)
    {
        return rate * (xold - xnew);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Exponential  Distribution : class  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Exponential  Distribution derived class. */
class Exponential_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Exponential_Distribution() : distr_gen(1.)
    {
        _init_defaults_("exp", {"rate"}, {1.},
                        0., NULL_DOUBLE, false, false);
    }
    ~Exponential_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));

        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0.) || (!tools::is_na(get_max())))
        {
            init_range(false, false);
            _init_truncation_(((discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE));
        }
    }

protected:
    std::exponential_distribution<double> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::exponential_distribution<double>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return exponential_impl::pdf(x, get_arg(0));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return exponential_impl::cdf(x, get_arg(0));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return exponential_impl::logpdf(x, get_arg(0));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return exponential_impl::logratio(xold, xnew, get_arg(0));
    }
    double _rand_(rnd_gen &prng) override
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
/////////////////////////////////  Poisson Distribution : methods  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace poisson_impl
{
    inline double pdf(double x, double lambda)
    {
        return exp(x * log(lambda) - lambda - lgamma(x + 1.));
    }
    inline double cdf(double x, double lambda)
    {
        return 1. - gamma_impl_core::pgamma(x + 1., lambda);
    }
    inline double logpdf(double x, double lambda)
    {
        return x * log(lambda) - lambda - lgamma(x + 1.);
    }
    inline double logratio(double xold, double xnew, double lambda)
    {
        return (xnew - xold) * log(lambda) + lgamma(xold + 1.) - lgamma(xnew + 1.);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Poisson Distribution : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson Distribution derived class. */
class Poisson_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Poisson_Distribution() : distr_gen(1.)
    {
        _init_defaults_("pois", {"lambda"}, {1.},
                        0., NULL_DOUBLE, false, false);
    }
    ~Poisson_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(D_DISCR_TYPE_DFND);

        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0.) || (!tools::is_na(get_max())))
        {
            init_range(false, false);
            _init_truncation_(D_DISCR_TYPE_DFND);
        }
    }

protected:
    std::poisson_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::poisson_distribution<unsigned>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return poisson_impl::pdf(x, get_arg(0));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return poisson_impl::cdf(x, get_arg(0));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return poisson_impl::logpdf(x, get_arg(0));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return poisson_impl::logratio(xold, xnew, get_arg(0));
    }
    double _rand_(rnd_gen &prng) override
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
////////////////////////////////  Binomial Distribution : methods  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace binomial_impl
{
    inline double pdf(double x, double size, double prob)
    {
        return exp(bax::log_binomial_coefficient(size, x) + x * log(prob) + (size - x) * log(1. - prob));
    }
    inline double cdf(double x, double size, double prob)
    {
        return beta_impl_core::ibeta(1. - prob, size - x, x + 1);
    }
    inline double logpdf(double x, double size, double prob)
    {
        return bax::log_binomial_coefficient(size, x) + x * log(prob) + (size - x) * log(1. - prob);
    }
    inline double logratio(double xold, double xnew, double size, double prob)
    {
        return bax::logratio_binomial_coefficient(size, xold, xnew) + (xnew - xold) * log(prob) + (xold - xnew) * log(1. - prob);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Binomial Distribution : class  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Binomial Distribution derived class. */
class Binomial_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Binomial_Distribution() : distr_gen(4, .5)
    {
        _init_defaults_("binom", {"size", "prob"}, {4, .5},
                        0., NULL_DOUBLE, false, false);
    }
    ~Binomial_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(D_DISCR_TYPE_DFND);

        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0.) || (!tools::is_na(get_max())))
        {
            init_range(false, false);
            _init_truncation_(D_DISCR_TYPE_DFND);
        }
    }

protected:
    std::binomial_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::binomial_distribution<unsigned>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return binomial_impl::pdf(x, get_arg(0), get_arg(1));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return binomial_impl::cdf(x, get_arg(0), get_arg(1));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return binomial_impl::logpdf(x, get_arg(0), get_arg(1));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return binomial_impl::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double _rand_(rnd_gen &prng) override
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
////////////////////////////  NegativeBinomial Distribution : methods  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace negativebinomial_impl
{
    inline double pdf(double x, double size, double prob)
    {
        return exp(bax::log_binomial_coefficient(x + size - 1, x) + x * log(1. - prob) + size * log(prob));
    }
    inline double cdf(double x, double size, double prob)
    {
        return beta_impl_core::ibeta(prob, size, x + 1);
    }
    inline double logpdf(double x, double size, double prob)
    {
        return bax::log_binomial_coefficient(x + size - 1, x) + x * log(1. - prob) + size * log(prob);
    }
    inline double logratio(double xold, double xnew, double size, double prob)
    {
        double output = 0.;
        // ratio of the binomial coefficients
        for (double i = 1.; i < size; ++i)
            output += log((i + xnew) / (i + xold));
        return output + (xnew - xold) * log(1. - prob);
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  NegativeBinomial Distribution : class  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* NegativeBinomial Distribution derived class. */
class NegativeBinomial_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    NegativeBinomial_Distribution() : distr_gen(4, .5)
    {
        _init_defaults_("nbinom", {"size", "prob"}, {4, .5},
                        0., NULL_DOUBLE, false, false);
    }
    ~NegativeBinomial_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(D_DISCR_TYPE_DFND);

        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0.) || (!tools::is_na(get_max())))
        {
            init_range(false, false);
            _init_truncation_(D_DISCR_TYPE_DFND);
        }
    }

protected:
    std::negative_binomial_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::negative_binomial_distribution<unsigned>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return binomial_impl::pdf(x, get_arg(0), get_arg(1));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return negativebinomial_impl::cdf(x, get_arg(0), get_arg(1));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return negativebinomial_impl::logpdf(x, get_arg(0), get_arg(1));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return negativebinomial_impl::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double _rand_(rnd_gen &prng) override
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
///////////////////////////////  Geometric Distribution : methods  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace geometric_impl
{
    inline double pdf(double x, double proba)
    {
        return proba * pow(1. - proba, x - 1.);
    }
    inline double cdf(double x, double proba)
    {
        return 1. - x * pow(1. - proba, x);
    }
    inline double logpdf(double x, double proba)
    {
        return log(proba) + (x - 1.) * log(1. - proba);
    }
    inline double logratio(double xold, double xnew, double proba)
    {
        return (xnew - xold) * log(1. - proba);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Geometric  Distribution : class  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Geometric Distribution derived class. */
class Geometric_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Geometric_Distribution() : distr_gen(1.)
    {
        _init_defaults_("geom", {"proba"}, {0.5},
                        0., NULL_DOUBLE, false, false);
    }
    ~Geometric_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        _init_discretization_(D_DISCR_TYPE_DFND);

        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0.) || (!tools::is_na(get_max())))
        {
            init_range(false, false);
            _init_truncation_(D_DISCR_TYPE_DFND);
        }
    }

protected:
    std::geometric_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void _update_args_() override
    {
        distr_gen.param(std::geometric_distribution<unsigned>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double _pdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return geometric_impl::pdf(x, get_arg(0));
    }
    double _cdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return geometric_impl::cdf(x, get_arg(0));
    }
    double _logpdf_(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return geometric_impl::logpdf(x, get_arg(0));
    }
    double _logratio_(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return geometric_impl::logratio(xold, xnew, get_arg(0));
    }
    double _rand_(rnd_gen &prng) override
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
        "unif", "norm", "lnorm", "gamma", "weibull", "rayleigh", "beta", "exp",
        "pois", "binom", "nbinom", "geom"};

    if (tools::is_na(name))
        return std::make_shared<Distribution>();

    else if (name == "unif")
        return std::make_shared<Uniform_Distribution>();

    else if (name == "norm")
        return std::make_shared<Normal_Distribution>();

    else if (name == "lnorm")
        return std::make_shared<Lognormal_Distribution>();

    else if (name == "gamma")
        return std::make_shared<Gamma_Distribution>();

    else if (name == "weibull")
        return std::make_shared<Weibull_Distribution>();

    else if (name == "rayleigh")
        return std::make_shared<Rayleigh_Distribution>();

    else if (name == "beta")
        return std::make_shared<Beta_Distribution>();

    else if (name == "exp")
        return std::make_shared<Exponential_Distribution>();

    else if (name == "pois")
        return std::make_shared<Poisson_Distribution>();

    else if (name == "binom")
        return std::make_shared<Binomial_Distribution>();

    else if (name == "nbinom")
        return std::make_shared<NegativeBinomial_Distribution>();

    else if (name == "geom")
        return std::make_shared<Geometric_Distribution>();

    else
        throw std::invalid_argument(
            "Unknown distribution \"" + name + "\"." + tools::find_closest(name, options));
}
