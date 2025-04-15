/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * distribution_base.cpp
 *
 */

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
    return (this->*_func_pdf)(x, opts);
}
double Distribution::cdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_cdf)(x, opts);
}
double Distribution::logpdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_logpdf)(x, opts);
}
double Distribution::logratio(double xold, double xnew, const dbl_vec &opts) const
{
    if (xold == xnew)
        return 0.;
    else
        return (this->*_func_logratio)(xold, xnew, opts);
}
double Distribution::rand(rnd_gen &prng)
{
    update_args();
    return (this->*_func_rand)(prng);
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
double Distribution::_rejection_sampling(rnd_gen &prng,
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

double Distribution::__default_pdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc)(raw_pdf(x, opts), opts);
}
double Distribution::__default_cdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc)((this->*_func_shift)(raw_cdf(x, opts), opts), opts);
}
double Distribution::__default_logpdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_logtrunc)(raw_logpdf(x, opts), opts);
}
double Distribution::__default_logratio(double xold, double xnew, const dbl_vec &opts) const
{
    return raw_logratio(xold, xnew, opts);
}
double Distribution::__default_rand(rnd_gen &prng)
{
    return raw_rand(prng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized defined ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__discr_defined_pdf(double x, const dbl_vec &opts) const
{
    return __default_pdf(std::floor(x), opts);
}
double Distribution::__discr_defined_cdf(double x, const dbl_vec &opts) const
{
    return __default_cdf(std::floor(x), opts);
}
double Distribution::__discr_defined_logpdf(double x, const dbl_vec &opts) const
{
    return __default_logpdf(std::floor(x), opts);
}
double Distribution::__discr_defined_logratio(double xold, double xnew, const dbl_vec &opts) const
{
    return __default_logratio(std::floor(xold), std::floor(xnew), opts);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized CENTER /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__discr_center_pdf(double x, const dbl_vec &opts) const
{
    const double d = round(x);
    return (this->*_func_trunc)(raw_cdf(d + .5, opts) - raw_cdf(d - .5, opts), opts);
}

double Distribution::__discr_center_cdf(double x, const dbl_vec &opts) const
{
    const double d = round(x);
    return (this->*_func_trunc)((this->*_func_shift)(raw_cdf(d + .5, opts), opts), opts);
}

double Distribution::__discr_center_logpdf(double x, const dbl_vec &opts) const
{
    return std::log(__discr_center_pdf(x, opts));
}

double Distribution::__discr_center_logratio(double xold, double xnew, const dbl_vec &opts) const
{
    const double dold = std::round(xold), dnew = std::round(xnew);
    return std::log((raw_cdf(dnew + .5, opts) - raw_cdf(dnew - .5, opts)) / (raw_cdf(dold + .5, opts) - raw_cdf(dold - .5, opts)));
}

double Distribution::__discr_center_rand(rnd_gen &prng) { return static_cast<double>(std::lround(__default_rand(prng))); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : discretized LEFT  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__discr_left_pdf(double x, const dbl_vec &opts) const
{
    const double d = std::floor(x);
    return (this->*_func_trunc)(raw_cdf(d + 1., opts) - raw_cdf(d, opts), opts);
}

double Distribution::__discr_left_cdf(double x, const dbl_vec &opts) const
{
    return (this->*_func_trunc)((this->*_func_shift)(raw_cdf(std::floor(x) + 1., opts), opts), opts);
}

double Distribution::__discr_left_logpdf(double x, const dbl_vec &opts) const
{
    return std::log(__discr_left_pdf(x, opts));
}

double Distribution::__discr_left_logratio(double xold, double xnew, const dbl_vec &opts) const
{
    const double dold = std::floor(xold), dnew = std::floor(xnew);
    return std::log((raw_cdf(dnew + 1., opts) - raw_cdf(dnew, opts)) / (raw_cdf(dold + 1., opts) - raw_cdf(dold, opts)));
}

double Distribution::__discr_left_rand(rnd_gen &prng)
{
    return std::floor(__default_rand(prng));
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  Base Distribution : truncation (default)  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__trunc_dummy(double x, const dbl_vec &opts) const
{
    UNUSED(opts);
    return x;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Base Distribution : truncation (continuos)  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__trunc_min(double x, const dbl_vec &opts) const
{
    return x / (1. - raw_cdf(get_min(), opts));
}
double Distribution::__trunc_max(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(get_max(), opts));
}
double Distribution::__trunc_both_sides(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(get_max(), opts) - raw_cdf(get_min(), opts));
}
double Distribution::__trunc_rand(rnd_gen &prng)
{
    return _rejection_sampling(prng, &Distribution::__default_rand);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  Base Distribution : truncation (discrete defined)  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__trunc_min_discr_defined(double x, const dbl_vec &opts) const
{
    return x / (1. - raw_cdf(get_min() - 1, opts));
}
double Distribution::__trunc_max_discr_defined(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(get_max() + 1, opts));
}
double Distribution::__trunc_both_sides_discr_defined(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(get_max() + 1, opts) - raw_cdf(get_min() - 1, opts));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  Base Distribution : truncation (CENTER discretization)  ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//----- Case of CENTER Discretization -----//
double Distribution::__trunc_min_discr_center(double x, const dbl_vec &opts) const
{
    return x / (1. - raw_cdf(std::ceil(get_min()) - .5, opts));
}
double Distribution::__trunc_max_discr_center(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(std::floor(get_max()) + .5, opts));
}
double Distribution::__trunc_both_sides_discr_center(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(std::floor(get_max()) + .5, opts) - raw_cdf(std::ceil(get_min()) - .5, opts));
}
double Distribution::__trunc_rand_discr_center(rnd_gen &prng)
{
    return _rejection_sampling(prng, &Distribution::__discr_center_rand);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Base Distribution : truncation (LEFT discretization)  /////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__trunc_min_discr_left(double x, const dbl_vec &opts) const
{
    return x / (1. - raw_cdf(std::floor(get_min()), opts));
}
double Distribution::__trunc_max_discr_left(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(std::ceil(get_max() + 1), opts));
}
double Distribution::__trunc_both_sides_discr_left(double x, const dbl_vec &opts) const
{
    return x / (raw_cdf(std::ceil(get_max() + 1), opts) - raw_cdf(std::floor(get_min()), opts));
}
double Distribution::__trunc_rand_discr_left(rnd_gen &prng)
{
    return _rejection_sampling(prng, &Distribution::__discr_left_rand);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Base Distribution : log-truncation (default)  /////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_dummy(double x, const dbl_vec &opts) const
{
    UNUSED(opts);
    return x;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Base Distribution : log-truncation (continuos)  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_min(double x, const dbl_vec &opts) const
{
    return x - std::log(1. - raw_cdf(get_min(), opts));
}
double Distribution::__logtrunc_max(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(get_max(), opts));
}
double Distribution::__logtrunc_both_sides(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(get_max(), opts) - raw_cdf(get_min(), opts));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Base Distribution : log-truncation (discrete defined)  ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_min_discr_defined(double x, const dbl_vec &opts) const
{
    return x - std::log(1. - raw_cdf(get_min() - 1, opts));
}
double Distribution::__logtrunc_max_discr_defined(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(get_max() + 1, opts));
}
double Distribution::__logtrunc_both_sides_discr_defined(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(get_max() + 1, opts) - raw_cdf(get_min() - 1, opts));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////  Base Distribution : log-truncation (CENTER discretization)  /////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_min_discr_center(double x, const dbl_vec &opts) const
{
    return x - std::log(1. - raw_cdf(std::ceil(get_min()) - .5, opts));
}
double Distribution::__logtrunc_max_discr_center(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(std::floor(get_max()) + .5, opts));
}
double Distribution::__logtrunc_both_sides_discr_center(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(std::floor(get_max()) + .5, opts) - raw_cdf(std::ceil(get_min()) - .5, opts));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  Base Distribution : log-truncation (LEFT discretization)  //////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_min_discr_left(double x, const dbl_vec &opts) const
{
    return x - std::log(1. - raw_cdf(std::floor(get_min()), opts));
}
double Distribution::__logtrunc_max_discr_left(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(std::ceil(get_max()), opts));
}
double Distribution::__logtrunc_both_sides_discr_left(double x, const dbl_vec &opts) const
{
    return x - std::log(raw_cdf(std::ceil(get_max()), opts) - raw_cdf(std::floor(get_min()), opts));
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  Base Distribution : shift cumulative (functional)  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__shift_dummy(double x, const dbl_vec &opts) const
{
    UNUSED(opts);
    return x;
}
double Distribution::__shift_min(double x, const dbl_vec &opts) const
{
    return x - raw_cdf(get_min(), opts);
}
double Distribution::__shift_min_discr_defined(double x, const dbl_vec &opts) const
{
    return x - raw_cdf(get_min() - 1, opts);
}
double Distribution::__shift_min_discr_center(double x, const dbl_vec &opts) const
{
    return x - raw_cdf(std::floor(get_min()) - .5, opts);
}
double Distribution::__shift_min_discr_left(double x, const dbl_vec &opts) const
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
///////////////////////////  Base Distribution : functional reshaping  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Tell if the distribution is discrete. */
bool Distribution::is_discrete() const { return _func_pdf != &Distribution::__default_pdf; }

/* Select functional method in case of discretization. */
void Distribution::init_functional_discretization(ushort discr_type)
{
    //===== DISCRETIZATION =====//
    switch (discr_type)
    {
    case D_DISCR_TYPE_DFND: // discrete defined
        _func_pdf = &Distribution::__discr_defined_pdf;
        _func_cdf = &Distribution::__discr_defined_cdf;
        _func_logpdf = &Distribution::__discr_defined_logpdf;
        _func_logratio = &Distribution::__discr_defined_logratio;
        _func_rand = &Distribution::__default_rand;
        break;

    case D_DISCR_TYPE_CNTR: // centered discretization
        _func_pdf = &Distribution::__discr_center_pdf;
        _func_cdf = &Distribution::__discr_center_cdf;
        _func_logpdf = &Distribution::__discr_center_logpdf;
        _func_logratio = &Distribution::__discr_center_logratio;
        _func_rand = &Distribution::__discr_center_rand;
        break;

    case D_DISCR_TYPE_LEFT: // left discretization
        _func_pdf = &Distribution::__discr_left_pdf;
        _func_cdf = &Distribution::__discr_left_cdf;
        _func_logpdf = &Distribution::__discr_left_logpdf;
        _func_logratio = &Distribution::__discr_left_logratio;
        _func_rand = &Distribution::__discr_left_rand;
        break;

    default: // no discretization (default)
        _func_pdf = &Distribution::__default_pdf;
        _func_cdf = &Distribution::__default_cdf;
        _func_logpdf = &Distribution::__default_logpdf;
        _func_logratio = &Distribution::__default_logratio;
        _func_rand = &Distribution::__default_rand;
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
        case D_DISCR_TYPE_DFND: // discrete defined
            _func_trunc = &Distribution::__trunc_min_discr_defined;
            _func_logtrunc = &Distribution::__logtrunc_min_discr_defined;
            _func_shift = &Distribution::__shift_min_discr_defined;
            _func_rand = &Distribution::__trunc_rand;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            _func_trunc = &Distribution::__trunc_min_discr_center;
            _func_logtrunc = &Distribution::__logtrunc_min_discr_center;
            _func_shift = &Distribution::__shift_min_discr_center;
            _func_rand = &Distribution::__trunc_rand_discr_center;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc = &Distribution::__trunc_min_discr_left;
            _func_logtrunc = &Distribution::__logtrunc_min_discr_left;
            _func_shift = &Distribution::__shift_min_discr_left;
            _func_rand = &Distribution::__trunc_rand_discr_left;
            break;
        default: // no discretization
            _func_trunc = &Distribution::__trunc_min;
            _func_logtrunc = &Distribution::__logtrunc_min;
            _func_shift = &Distribution::__shift_min;
            _func_rand = &Distribution::__trunc_rand;
        }
        break;

    case 2: //----- upper range truncation -----//
        _func_shift = &Distribution::__shift_dummy;
        switch (discr_type)
        {
        case D_DISCR_TYPE_DFND: // discrete defined
            _func_trunc = &Distribution::__trunc_max_discr_defined;
            _func_logtrunc = &Distribution::__logtrunc_max_discr_defined;
            _func_rand = &Distribution::__trunc_rand;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            _func_trunc = &Distribution::__trunc_max_discr_center;
            _func_logtrunc = &Distribution::__logtrunc_max_discr_center;
            _func_rand = &Distribution::__trunc_rand_discr_center;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc = &Distribution::__trunc_max_discr_left;
            _func_logtrunc = &Distribution::__logtrunc_max_discr_left;
            _func_rand = &Distribution::__trunc_rand_discr_left;
            break;
        default: // no discretization
            _func_trunc = &Distribution::__trunc_max;
            _func_logtrunc = &Distribution::__logtrunc_max;
            _func_rand = &Distribution::__trunc_rand;
        }
        break;

    case 3: //----- both sides truncation -----//
        switch (discr_type)
        {
        case D_DISCR_TYPE_DFND: // discrete defined
            _func_trunc = &Distribution::__trunc_both_sides_discr_defined;
            _func_logtrunc = &Distribution::__logtrunc_both_sides_discr_defined;
            _func_shift = &Distribution::__shift_min_discr_defined;
            _func_rand = &Distribution::__trunc_rand;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            _func_trunc = &Distribution::__trunc_both_sides_discr_center;
            _func_logtrunc = &Distribution::__logtrunc_both_sides_discr_center;
            _func_shift = &Distribution::__shift_min_discr_center;
            _func_rand = &Distribution::__trunc_rand_discr_center;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _func_trunc = &Distribution::__trunc_both_sides_discr_left;
            _func_logtrunc = &Distribution::__logtrunc_both_sides_discr_left;
            _func_shift = &Distribution::__shift_min_discr_left;
            _func_rand = &Distribution::__trunc_rand_discr_left;
            break;
        default: // no discretization
            _func_trunc = &Distribution::__trunc_both_sides;
            _func_logtrunc = &Distribution::__logtrunc_both_sides;
            _func_shift = &Distribution::__shift_min;
            _func_rand = &Distribution::__trunc_rand;
        }
        break;

    default: //----- no truncation (case 0) -----//
        _func_trunc = &Distribution::__trunc_dummy;
        _func_logtrunc = &Distribution::__logtrunc_dummy;
        _func_shift = &Distribution::__shift_dummy;
        break;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Base Distribution - reshaping  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Default reshaping initialization.
 * Necessary to have dynamic state variables for NULL case
 */
void Distribution::init_reshaping(bool discrete)
{
    ushort discr_case = (discrete) ? D_DISCR_TYPE_DFND : D_DISCR_TYPE_NONE;
    init_functional_discretization(discr_case);
    init_functional_truncation(discr_case);
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
    init_functional_discretization(D_DISCR_TYPE_DFND);

    if (get_min() < 0.)
        throw std::invalid_argument(
            "Invalid min value " + show_min() + " for distribution \"" + show_name() + "\".");
    else if ((get_min() > 0.) || (!tools::is_na(get_max())))
    {
        init_functional_range(false, false);
        init_functional_truncation(D_DISCR_TYPE_DFND);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Distribution : showcase (method)  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Return values to showcase the profile of the distribution. */
const std::vector<dbl_vec> Distribution::showcase(const dbl_vec &opts,
                                                  int n_bins,
                                                  double ytol) const
{
    if (n_bins < 1)
        throw std::invalid_argument("Number of bins should be greater than 0.");

    // pick extremes for the showcase
    double lower = get_min();
    double upper = get_max();

    if (tools::is_na(upper))
    {
        upper = 10. + (tools::is_na(lower) ? 0 : lower);
        while (1. - cdf(upper, opts) > ytol)
            ++upper;
    }

    if (tools::is_na(lower))
    {
        lower = -10.;
        while (cdf(lower, opts) > ytol)
            --lower;
    }

    // compute binning width
    double delta_x = is_discrete() ? 1. : ((upper - lower) / n_bins);

    // adjust extremes if interval is open
    if (!valid(lower))
        lower += delta_x;
    if (!valid(upper))
        upper -= delta_x;

    // re-compute binning width
    delta_x = is_discrete() ? 1. : ((upper - lower) / n_bins);

    size_t N = (upper - lower) / delta_x;
    dbl_vec gentime(N + 1);
    dbl_vec inst(N + 1);
    dbl_vec cumul(N + 1);

    for (size_t i = 0; i <= N; ++i)
    {
        double x = lower + i * delta_x;
        gentime[i] = x;
        inst[i] = pdf(x, opts);
        cumul[i] = cdf(x, opts);
    }

    return {gentime, inst, cumul};
}

/* Developer method to check for normalization. */
double Distribution::__d_e_v__check__pdf__norm(const dbl_vec &opts,
                                                  int n_bins,
                                                  double ytol) const
{
    // produce a showcase of the distribution profile
    const std::vector<dbl_vec> values = showcase(opts, n_bins, ytol);

    // compute mean and variance
    const dbl_vec &x = values[0];
    const dbl_vec &prob = values[1];
    const size_t N = x.size();
    const double Dx = x[1] - x[0];

    // note: integration with trapezoidal rule for continuos case
    double m = prob[0] + prob[N - 1];
    if (is_discrete())
        m *= 0.5;
    for (size_t i = 1; i < N - 1; ++i)
        m += prob[i];
    m *= Dx;

    return m;
}
