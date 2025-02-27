/*
 * distribution.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "numerical_recipes.h"
#include "distribution.h"

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
    if (xold == xnew)
        return 0.;
    else
        return (this->*_func_logratio_)(xold, xnew, opts);
}
double Distribution::rand(rnd_gen &prng)
{
    update_args();
    return (this->*_func_rand_)(prng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution - virtual  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::raw_pdf(double x, const dbl_vec &opts) const
{
    UNUSED(x);
    UNUSED(opts);
    throw std::invalid_argument(
        "Attempt to call undefined method \"pdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::raw_cdf(double x, const dbl_vec &opts) const
{
    UNUSED(x);
    UNUSED(opts);
    throw std::invalid_argument(
        "Attempt to call undefined method \"cdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::raw_logpdf(double x, const dbl_vec &opts) const
{
    UNUSED(x);
    UNUSED(opts);
    throw std::invalid_argument(
        "Attempt to call undefined method \"logpdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::raw_logratio(double xold, double xnew, const dbl_vec &opts) const
{
    UNUSED(xold);
    UNUSED(xnew);
    UNUSED(opts);
    return 0.;
}
void Distribution::update_args() {};
double Distribution::raw_rand(rnd_gen &prng)
{
    UNUSED(prng);
    throw std::invalid_argument(
        "Attempt to call undefined method \"rand\" of distribution \"" + show_name() + "\".");
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Rejection sampling. */
double Distribution::_rejection_sampling_(rnd_gen &prng,
                                          double (Distribution::*func)(rnd_gen &))
{
    double xnew = (this->*func)(prng);

    unsigned n_attempts = 0;
    while (!valid(xnew))
    {
        if (n_attempts > MAX_NUM_ATTEMPTS)
            throw std::runtime_error(
                "Exceeded attempts (N > " + std::to_string(MAX_NUM_ATTEMPTS) + ") in random draw of " + show_name() + "(" + show_arg_values() + ")" + " with truncation " + show_range() + ".");
        xnew = (this->*func)(prng);
        ++n_attempts;
    }

    return xnew;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : defaults  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_default_pdf_(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc_)(raw_pdf(x, opts), opts);
}
double Distribution::_default_cdf_(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc_)((this->*_func_shift_)(raw_cdf(x, opts), opts), opts);
}
double Distribution::_default_logpdf_(double x, const dbl_vec &opts) const
{
    return (this->*_func_logtrunc_)(raw_logpdf(x, opts), opts);
}
double Distribution::_default_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    return raw_logratio(xold, xnew, opts);
}
double Distribution::_default_rand_(rnd_gen &prng)
{
    return raw_rand(prng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized defined ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_discr_defined_pdf_(double x, const dbl_vec &opts) const
{
    return _default_pdf_(std::floor(x), opts);
}
double Distribution::_discr_defined_cdf_(double x, const dbl_vec &opts) const
{
    return _default_cdf_(std::floor(x), opts);
}
double Distribution::_discr_defined_logpdf_(double x, const dbl_vec &opts) const
{
    return _default_logpdf_(std::floor(x), opts);
}
double Distribution::_discr_defined_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    return _default_logratio_(std::floor(xold), std::floor(xnew), opts);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized CENTER /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_discr_center_pdf_(double x, const dbl_vec &opts) const
{
    const double d = round(x);
    return (this->*_func_trunc_)(raw_cdf(d + .5, opts) - raw_cdf(d - .5, opts), opts);
}

double Distribution::_discr_center_cdf_(double x, const dbl_vec &opts) const
{
    const double d = round(x);
    return (this->*_func_trunc_)((this->*_func_shift_)(raw_cdf(d + .5, opts), opts), opts);
}

double Distribution::_discr_center_logpdf_(double x, const dbl_vec &opts) const
{
    return std::log(_discr_center_pdf_(x, opts));
}

double Distribution::_discr_center_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    const double dold = std::round(xold), dnew = std::round(xnew);
    return std::log((raw_cdf(dnew + .5, opts) - raw_cdf(dnew - .5, opts)) / (raw_cdf(dold + .5, opts) - raw_cdf(dold - .5, opts)));
}

double Distribution::_discr_center_rand_(rnd_gen &prng)
{
    // UPDATE : lround to avoid casting like round(-0.1) -> -0.0
    return static_cast<double>(std::lround(_default_rand_(prng)));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : discretized LEFT  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_discr_left_pdf_(double x, const dbl_vec &opts) const
{
    const double d = std::floor(x);
    return (this->*_func_trunc_)(raw_cdf(d + 1., opts) - raw_cdf(d, opts), opts);
}

double Distribution::_discr_left_cdf_(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc_)((this->*_func_shift_)(raw_cdf(std::floor(x) + 1., opts), opts), opts);
}

double Distribution::_discr_left_logpdf_(double x, const dbl_vec &opts) const
{
    return std::log(_discr_left_pdf_(x, opts));
}

double Distribution::_discr_left_logratio_(double xold, double xnew, const dbl_vec &opts) const
{
    const double dold = std::floor(xold), dnew = std::floor(xnew);
    return std::log((raw_cdf(dnew + 1., opts) - raw_cdf(dnew, opts)) / (raw_cdf(dold + 1., opts) - raw_cdf(dold, opts)));
}

double Distribution::_discr_left_rand_(rnd_gen &prng)
{
    return std::floor(_default_rand_(prng));
}

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
    return x / (1. - raw_cdf(get_min(), opts));
}
double Distribution::_trunc_max_(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(get_max(), opts));
}
double Distribution::_trunc_both_sides_(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(get_max(), opts) - raw_cdf(get_min(), opts));
}
double Distribution::_trunc_rand_(rnd_gen &prng)
{
    return _rejection_sampling_(prng, &Distribution::_default_rand_);
}

//----- Case of CENTER Discretization -----//
double Distribution::_trunc_min_discr_center_(double x, const dbl_vec &opts) const
{
    return x / (1. - raw_cdf(std::ceil(get_min()) - .5, opts));
}
double Distribution::_trunc_max_discr_center_(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(std::floor(get_max()) + .5, opts));
}
double Distribution::_trunc_both_sides_discr_center_(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(std::floor(get_max()) + .5, opts) - raw_cdf(std::ceil(get_min()) - .5, opts));
}
double Distribution::_trunc_rand_discr_center_(rnd_gen &prng)
{
    return _rejection_sampling_(prng, &Distribution::_discr_center_rand_);
}

//----- Case of LEFT Discretization -----//
double Distribution::_trunc_min_discr_left_(double x, const dbl_vec &opts) const
{
    return x / (1. - raw_cdf(std::floor(get_min()), opts));
}
double Distribution::_trunc_max_discr_left_(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(std::ceil(get_max()), opts));
}
double Distribution::_trunc_both_sides_discr_left_(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(std::ceil(get_max()), opts) - raw_cdf(std::floor(get_min()), opts));
}
double Distribution::_trunc_rand_discr_left_(rnd_gen &prng)
{
    return _rejection_sampling_(prng, &Distribution::_discr_left_rand_);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Base Distribution : logtruncation (functional)  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_logtrunc_dummy_(double x, const dbl_vec &opts) const
{
    UNUSED(opts);
    return x;
}
double Distribution::_logtrunc_min_(double x, const dbl_vec &opts) const
{
    return x - std::log(1. - raw_cdf(get_min(), opts));
}
double Distribution::_logtrunc_max_(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(get_max(), opts));
}
double Distribution::_logtrunc_both_sides_(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(get_max(), opts) - raw_cdf(get_min(), opts));
}

//----- Case of Center Discretization -----//
double Distribution::_logtrunc_min_discr_center_(double x, const dbl_vec &opts) const
{
    return x - std::log(1. - raw_cdf(std::ceil(get_min()) - .5, opts));
}
double Distribution::_logtrunc_max_discr_center_(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(std::floor(get_max()) + .5, opts));
}
double Distribution::_logtrunc_both_sides_discr_center_(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(std::floor(get_max()) + .5, opts) - raw_cdf(std::ceil(get_min()) - .5, opts));
}

//----- Case of Left Discretization -----//
double Distribution::_logtrunc_min_discr_left_(double x, const dbl_vec &opts) const
{
    return x - std::log(1. - raw_cdf(std::floor(get_min()), opts));
}
double Distribution::_logtrunc_max_discr_left_(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(std::ceil(get_max()), opts));
}
double Distribution::_logtrunc_both_sides_discr_left_(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(std::ceil(get_max()), opts) - raw_cdf(std::floor(get_min()), opts));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  Base Distribution : shift cumulative (functional)  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::_shift_dummy_(double x, const dbl_vec &opts) const
{
    UNUSED(opts);
    return x;
}
double Distribution::_shift_min_(double x, const dbl_vec &opts) const
{
    return x - raw_cdf(get_min(), opts);
}
double Distribution::_shift_min_discr_center_(double x, const dbl_vec &opts) const
{
    return x - raw_cdf(std::floor(get_min()) - .5, opts);
}
double Distribution::_shift_min_discr_left_(double x, const dbl_vec &opts) const
{
    return x - raw_cdf(std::floor(get_min()), opts);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Base Distribution - reshaping  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Default reshaping initialization. */
void Distribution::init_reshaping(bool discrete)
{
    UNUSED(discrete);
}

/* Frequent reshaping while domain has zero lower bound. */
void Distribution::init_min_zero_reshaping(bool discrete, bool is_open_zero)
{
    ushort discr_case = (discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE;
    init_functional_discretization(discr_case);

    if (get_min() < 0.)
        throw std::invalid_argument(
            "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
    else if ((get_min() > 0.) || (!tools::is_na(get_max())))
    {
        init_functional_range(((is_open_zero) ? (get_min() == 0.) : false), false);
        init_functional_truncation(discr_case);
    }
}

/* Frequent reshaping while domain has zero lower bound for discrete defined distributions. */
void Distribution::init_min_zero_discr_defined_reshaping()
{
    if (get_min() < 0.)
        throw std::invalid_argument(
            "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
    else if ((get_min() > 0.) || (!tools::is_na(get_max())))
    {
        init_functional_range(false, false);
        init_functional_truncation(D_DISCR_TYPE_DFND);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : assign reshaping  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Tell if the distribution is discrete. */
bool Distribution::is_discrete() const { return _func_pdf_ != &Distribution::_default_pdf_; }

/* Select functional method in case of discretization. */
void Distribution::init_functional_discretization(ushort discr_type)
{
    //===== DISCRETIZATION =====//
    switch (discr_type)
    {
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

    default: // no discretization (default)
        _func_pdf_ = &Distribution::_default_pdf_;
        _func_cdf_ = &Distribution::_default_cdf_;
        _func_logpdf_ = &Distribution::_default_logpdf_;
        _func_logratio_ = &Distribution::_default_logratio_;
        _func_rand_ = &Distribution::_default_rand_;
    }
}

/* Select functional method in case of truncation. */
void Distribution::init_functional_truncation(ushort discr_type)
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
            _func_rand_ = &Distribution::_trunc_rand_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc_ = &Distribution::_trunc_min_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_min_discr_left_;
            _func_shift_ = &Distribution::_shift_min_discr_left_;
            _func_rand_ = &Distribution::_trunc_rand_discr_left_;
            break;
        default: // no discretization or D_DISCR_TYPE_DFND
            _func_trunc_ = &Distribution::_trunc_min_;
            _func_logtrunc_ = &Distribution::_logtrunc_min_;
            _func_shift_ = &Distribution::_shift_min_;
            _func_rand_ = &Distribution::_trunc_rand_;
        }
        break;

    case 2: //----- upper range truncation -----//
        _func_shift_ = &Distribution::_shift_dummy_;
        switch (discr_type)
        {
        case D_DISCR_TYPE_CNTR: // centered discretization
            _func_trunc_ = &Distribution::_trunc_max_discr_center_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_discr_center_;
            _func_rand_ = &Distribution::_trunc_rand_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc_ = &Distribution::_trunc_max_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_discr_left_;
            _func_rand_ = &Distribution::_trunc_rand_discr_left_;
            break;
        default: // no discretization or D_DISCR_TYPE_DFND
            _func_trunc_ = &Distribution::_trunc_max_;
            _func_logtrunc_ = &Distribution::_logtrunc_max_;
            _func_rand_ = &Distribution::_trunc_rand_;
        }
        break;

    case 3: //----- both sides truncation -----//
        switch (discr_type)
        {
        case D_DISCR_TYPE_CNTR: // centered discretization
            _func_trunc_ = &Distribution::_trunc_both_sides_discr_center_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_discr_center_;
            _func_shift_ = &Distribution::_shift_min_discr_center_;
            _func_rand_ = &Distribution::_trunc_rand_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc_ = &Distribution::_trunc_both_sides_discr_left_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_discr_left_;
            _func_shift_ = &Distribution::_shift_min_discr_left_;
            _func_rand_ = &Distribution::_trunc_rand_discr_left_;
            break;
        default: // no discretization or D_DISCR_TYPE_DFND
            _func_trunc_ = &Distribution::_trunc_both_sides_;
            _func_logtrunc_ = &Distribution::_logtrunc_both_sides_;
            _func_shift_ = &Distribution::_shift_min_;
            _func_rand_ = &Distribution::_trunc_rand_;
        }
        break;

    default: //----- no truncation (case 0) -----//
        _func_trunc_ = &Distribution::_trunc_dummy_;
        _func_logtrunc_ = &Distribution::_logtrunc_dummy_;
        _func_shift_ = &Distribution::_shift_dummy_;
        _func_rand_ = &Distribution::_trunc_rand_;
        break;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Distribution : showcase (method)  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return values to showcase the profile of the distribution. */
const std::map<double, dbl_vec> Distribution::showcase(const dbl_vec &opts,
                                                       const dbl_vec &xlim,
                                                       int n_bins) const
{
    if (n_bins < 1)
        throw std::invalid_argument("Number of bins should be greater than 0.");

    // pick extremes for the showcase
    double lower = get_min();
    double upper = get_max();

    // adjust lower bound and limit
    if (tools::is_na(lower))
        lower = xlim[0];
    else
        lower = (xlim[0] < lower) ? lower : xlim[0];

    // adjust upper bound and limit
    if (tools::is_na(upper))
        upper = (xlim[1] > lower) ? xlim[1] : lower + 20.;
    else
        upper = (xlim[1] > upper) ? upper : xlim[1];

    // adjust binning width
    double delta_x = is_discrete() ? 1. : ((upper - lower) / n_bins);

    // adjust extremes if interval is open
    if (!valid(lower))
        lower += delta_x;
    if (!valid(upper))
        upper -= delta_x;

    // compute values to showcase
    std::map<double, dbl_vec> output = {};
    for (double x = lower; x <= upper; x += delta_x)
        output.emplace(x, std::vector<double>({pdf(x, opts), cdf(x, opts)}));

    return output;
}

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
        init_defaults("unif", {"min", "max"}, {0., 1.},
                      NULL_DOUBLE, NULL_DOUBLE, false, false);
    }
    ~Uniform_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        ushort discr_case = (discrete) ? D_DISCR_TYPE_CNTR : D_DISCR_TYPE_NONE;

        init_functional_range(false, false);
        init_functional_discretization(discr_case);
        init_functional_truncation(discr_case);
    }

protected:
    std::uniform_real_distribution<double> distr_gen;

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(x);
        UNUSED(opts);
        return 1.;
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return x;
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(x);
        UNUSED(opts);
        return 0.;
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(xold);
        UNUSED(xnew);
        UNUSED(opts);
        return 0.;
    }
    double raw_rand(rnd_gen &prng) override
    {
        const double r = distr_gen(prng);
        return get_arg(0) * (1 - r) + get_arg(1) * r;
    }
};

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
        init_defaults("norm", {"mean", "sd"}, {0., 1.},
                      NULL_DOUBLE, NULL_DOUBLE, false, false);
    }
    ~Normal_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        ushort discr_case = (discrete) ? D_DISCR_TYPE_CNTR : D_DISCR_TYPE_NONE;

        init_functional_range(false, false);
        init_functional_discretization(discr_case);
        init_functional_truncation(discr_case);
    }

protected:
    std::normal_distribution<double> distr_gen;

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::normal::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::normal::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::normal::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::normal::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override
    {
        return get_arg(0) + get_arg(1) * distr_gen(prng);
    }
};

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
        init_defaults("lnorm", {"meanlog", "sdlog"}, {0., 1.},
                      0., NULL_DOUBLE, false, false);
    }
    ~Lognormal_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::lognormal_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override
    {
        distr_gen.param(std::lognormal_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::lognormal::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::lognormal::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::lognormal::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::lognormal::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override
    {
        return distr_gen(prng);
    }
};

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
        init_defaults("beta", {"shape1", "shape2"}, {1., 1.},
                      0., 1., true, true);
    }
    ~Beta_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        if (discrete)
            throw std::invalid_argument(
                "Invalid option discrete for distribution \"" + show_name() + "\".");
        init_functional_discretization(D_DISCR_TYPE_NONE);

        if (get_min() < 0)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if (get_max() > 1)
            throw std::invalid_argument(
                "Invalid max value " + show_max() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0) || (get_max() < 1))
        {
            init_functional_range((get_min() == 0), (get_max() == 1));
            init_functional_truncation(D_DISCR_TYPE_NONE);
        }
    }

protected:
    std::gamma_distribution<double> gamma_1;
    std::gamma_distribution<double> gamma_2;

    //==== SETTER =====//
    void update_args() override
    {
        gamma_1.param(std::gamma_distribution<double>::param_type(get_arg(0), 1.));
        gamma_2.param(std::gamma_distribution<double>::param_type(get_arg(1), 1.));
    }

    //==== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::beta::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::beta::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::beta::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::beta::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override
    {
        double x1 = gamma_1(prng);
        double x2 = gamma_2(prng);
        return x1 / (x1 + x2);
    }
};

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
        init_defaults("gamma", {"shape", "scale"}, {1., 1.},
                      0., NULL_DOUBLE, false, false);
    }
    ~Gamma_Distribution() = default;

    //==== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::gamma_distribution<double> distr_gen;

    //==== SETTER =====//
    void update_args() override
    {
        distr_gen.param(std::gamma_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //==== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::gamma::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::gamma::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::gamma::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::gamma::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override
    {
        return distr_gen(prng);
    }
};

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
        init_defaults("rayleigh", {"scale"}, {1.},
                      0., NULL_DOUBLE, false, false);
    }
    ~Rayleigh_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::uniform_real_distribution<double> unif_gen;

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::rayleigh::pdf(x, get_arg(0));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::rayleigh::cdf(x, get_arg(0));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::rayleigh::logpdf(x, get_arg(0));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::rayleigh::logratio(xold, xnew, get_arg(0));
    }
    double raw_rand(rnd_gen &prng) override
    {
        return get_arg(0) * sqrt(-2 * std::log(unif_gen(prng)));
    }
};

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
        init_defaults("weibull", {"shape", "scale"}, {1., 1.},
                      0., NULL_DOUBLE, false, false);
    }
    ~Weibull_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::weibull_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override
    {
        distr_gen.param(std::weibull_distribution<double>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::weibull::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::weibull::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::weibull::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::weibull::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override
    {
        return distr_gen(prng);
    }
};

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
        init_defaults("exp", {"rate"}, {1.},
                      0., NULL_DOUBLE, false, false);
    }
    ~Exponential_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        init_min_zero_reshaping(discrete, false);
    }

protected:
    std::exponential_distribution<double> distr_gen;

    //===== SETTER =====//
    void update_args() override
    {
        distr_gen.param(std::exponential_distribution<double>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::exponential::pdf(x, get_arg(0));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::exponential::cdf(x, get_arg(0));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::exponential::logpdf(x, get_arg(0));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::exponential::logratio(xold, xnew, get_arg(0));
    }
    double raw_rand(rnd_gen &prng) override
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
/////////////////////////////////  Poisson Distribution : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson Distribution derived class. */
class Poisson_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Poisson_Distribution() : distr_gen(1.)
    {
        init_defaults("pois", {"lambda"}, {1.},
                      0., NULL_DOUBLE, false, false);
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~Poisson_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        UNUSED(discrete);
        init_min_zero_discr_defined_reshaping();
    }

protected:
    std::poisson_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override
    {
        distr_gen.param(std::poisson_distribution<unsigned>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::poisson::pdf(x, get_arg(0));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::poisson::cdf(x, get_arg(0));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::poisson::logpdf(x, get_arg(0));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::poisson::logratio(xold, xnew, get_arg(0));
    }
    double raw_rand(rnd_gen &prng) override
    {
        return distr_gen(prng);
    }
};

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
        init_defaults("binom", {"size", "prob"}, {4, .5},
                      0., NULL_DOUBLE, false, false);
        init_functional_discretization(D_DISCR_TYPE_DFND);
        
        assign_pointer_range_max(get_ptr_to_arg(0));
        init_functional_range(false, false);
    }
    ~Binomial_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        UNUSED(discrete);

        if (get_min() < 0.)
            throw std::invalid_argument(
                "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
        else if ((get_min() > 0.) || (get_ptr_to_arg(0) != get_ptr_to_max()))
        {
            if (get_max() > get_arg(0))
                throw std::invalid_argument(
                    "Invalid max value " + show_max() + " for distribution \"" + show_name() + "\" with size " + std::to_string(get_arg(0)) + ".");
            init_functional_range(false, false);
            init_functional_truncation(D_DISCR_TYPE_DFND);
        }
    }

protected:
    std::binomial_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override
    {
        distr_gen.param(std::binomial_distribution<unsigned>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::binomial::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::binomial::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::binomial::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::binomial::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override
    {
        return distr_gen(prng);
    }
};

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
        init_defaults("nbinom", {"size", "prob"}, {4, .5},
                      0., NULL_DOUBLE, false, false);
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~NegativeBinomial_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        UNUSED(discrete);
        init_min_zero_discr_defined_reshaping();
    }

protected:
    std::negative_binomial_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override
    {
        distr_gen.param(std::negative_binomial_distribution<unsigned>::param_type(get_arg(0), get_arg(1)));
    }

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::negbinomial::pdf(x, get_arg(0), get_arg(1));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::negbinomial::cdf(x, get_arg(0), get_arg(1));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::negbinomial::logpdf(x, get_arg(0), get_arg(1));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::negbinomial::logratio(xold, xnew, get_arg(0), get_arg(1));
    }
    double raw_rand(rnd_gen &prng) override
    {
        return distr_gen(prng);
    }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Geometric Distribution : class  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Geometric Distribution derived class. */
class Geometric_Distribution : public Distribution
{
public:
    //===== CONSTRUCTOR =====//
    Geometric_Distribution() : distr_gen(1.)
    {
        init_defaults("geom", {"proba"}, {0.5},
                      0., NULL_DOUBLE, false, false);
        init_functional_discretization(D_DISCR_TYPE_DFND);
    }
    ~Geometric_Distribution() = default;

    //===== SETTER =====//
    void init_reshaping(bool discrete) override
    {
        UNUSED(discrete);
        init_min_zero_discr_defined_reshaping();
    }

protected:
    std::geometric_distribution<unsigned> distr_gen;

    //===== SETTER =====//
    void update_args() override
    {
        distr_gen.param(std::geometric_distribution<unsigned>::param_type(get_arg(0)));
    }

    //===== GETTER =====//
    double raw_pdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::geometric::pdf(x, get_arg(0));
    }
    double raw_cdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::geometric::cdf(x, get_arg(0));
    }
    double raw_logpdf(double x, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::geometric::logpdf(x, get_arg(0));
    }
    double raw_logratio(double xold, double xnew, const dbl_vec &opts) const override
    {
        UNUSED(opts);
        return numr::geometric::logratio(xold, xnew, get_arg(0));
    }
    double raw_rand(rnd_gen &prng) override
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
//////////////////////////////////  Distribution Switchboard  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Auxiliary hidden template to return the chosen distribution. */
std::shared_ptr<Distribution> distribution_switchboard(const std::string &name)
{
    const str_vec &options = {"unif", "norm", "lnorm",
                              "beta", "gamma", "rayleigh", "weibull", "exp",
                              "pois", "binom", "nbinom", "geom"};

    if (tools::is_na(name))
        return std::make_shared<Distribution>();

    else if (name == "unif")
        return std::make_shared<Uniform_Distribution>();

    else if (name == "norm")
        return std::make_shared<Normal_Distribution>();

    else if (name == "lnorm")
        return std::make_shared<Lognormal_Distribution>();

    else if (name == "beta")
        return std::make_shared<Beta_Distribution>();

    else if (name == "gamma")
        return std::make_shared<Gamma_Distribution>();

    else if (name == "rayleigh")
        return std::make_shared<Rayleigh_Distribution>();

    else if (name == "weibull")
        return std::make_shared<Weibull_Distribution>();

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
