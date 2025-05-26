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
//////////////////////////////  Base Distribution - constructor  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Distribution::Distribution(const std::string &name,
                           const str_vec &arg_name,
                           const dbl_vec &arg_value,
                           double domain_min,
                           double domain_max,
                           bool is_min_open,
                           bool is_max_open)
    // each distribution has extra argument "x" for return_value
    : Function(name, arg_name + str_vec({"x"}), arg_value + dbl_vec({HOUSEPI_NA_DOUBLE}),
               domain_min, domain_max, is_min_open, is_max_open)
{
}
Distribution::~Distribution() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution - wrapper  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::pdf(double x) const
{
    return (this->*_do_pdf)(x);
}
double Distribution::cdf(double x) const
{
    return (this->*_do_cdf)(x);
}
double Distribution::logpdf(double x) const
{
    return (this->*_do_logpdf)(x);
}
double Distribution::logratio(double xold, double xnew) const
{
    if (xold == xnew)
        return 0.;
    else
        return (this->*_do_logratio)(xold, xnew);
}
double Distribution::rand(rnd_gen &prng)
{
    update_args();
    return (this->*_do_rand)(prng);
}
double Distribution::return_value() const
{
    const double x = get_arg(n_args() - 1);
    if (valid(x))
        return pdf(x);
    else
        return 0.;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution - virtual  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::raw_pdf(double x) const
{
    UNUSED(x);
    throw std::invalid_argument(
        "Attempt to call undefined method \"pdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::raw_cdf(double x) const
{
    UNUSED(x);
    throw std::invalid_argument(
        "Attempt to call undefined method \"cdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::raw_logpdf(double x) const
{
    UNUSED(x);

    throw std::invalid_argument(
        "Attempt to call undefined method \"logpdf\" of distribution \"" + show_name() + "\".");
}
double Distribution::raw_logratio(double xold, double xnew) const
{
    UNUSED(xold);
    UNUSED(xnew);

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
                "Exceeded attempts (N > " + std::to_string(MAX_NUM_ATTEMPTS) + ") in random draw of " + show_summary() + ".");
        xnew = (this->*func)(prng);
        ++n_attempts;
    }

    return xnew;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : defaults  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__default_pdf(double x) const
{
    return (this->*_do_trunc)(raw_pdf(x));
}
double Distribution::__default_cdf(double x) const
{
    return (this->*_do_trunc)((this->*_do_shift)(raw_cdf(x)));
}
double Distribution::__default_logpdf(double x) const
{
    return (this->*_do_logtrunc)(raw_logpdf(x));
}
double Distribution::__default_logratio(double xold, double xnew) const
{
    return raw_logratio(xold, xnew);
}
double Distribution::__default_rand(rnd_gen &prng)
{
    return raw_rand(prng);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized defined ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__discr_defined_pdf(double x) const
{
    return __default_pdf(std::floor(x));
}
double Distribution::__discr_defined_cdf(double x) const
{
    return __default_cdf(std::floor(x));
}
double Distribution::__discr_defined_logpdf(double x) const
{
    return __default_logpdf(std::floor(x));
}
double Distribution::__discr_defined_logratio(double xold, double xnew) const
{
    return __default_logratio(std::floor(xold), std::floor(xnew));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Base Distribution : discretized CENTER /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__discr_center_pdf(double x) const
{
    const double d = round(x);
    return (this->*_do_trunc)(raw_cdf(d + .5) - raw_cdf(d - .5));
}

double Distribution::__discr_center_cdf(double x) const
{
    const double d = round(x);
    return (this->*_do_trunc)((this->*_do_shift)(raw_cdf(d + .5)));
}

double Distribution::__discr_center_logpdf(double x) const
{
    return std::log(__discr_center_pdf(x));
}

double Distribution::__discr_center_logratio(double xold, double xnew) const
{
    const double dold = std::round(xold), dnew = std::round(xnew);
    return std::log((raw_cdf(dnew + .5) - raw_cdf(dnew - .5)) / (raw_cdf(dold + .5) - raw_cdf(dold - .5)));
}

double Distribution::__discr_center_rand(rnd_gen &prng) { return static_cast<double>(std::lround(__default_rand(prng))); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : discretized LEFT  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__discr_left_pdf(double x) const
{
    const double d = std::floor(x);
    return (this->*_do_trunc)(raw_cdf(d + 1.) - raw_cdf(d));
}

double Distribution::__discr_left_cdf(double x) const
{
    return (this->*_do_trunc)((this->*_do_shift)(raw_cdf(std::floor(x) + 1.)));
}

double Distribution::__discr_left_logpdf(double x) const
{
    return std::log(__discr_left_pdf(x));
}

double Distribution::__discr_left_logratio(double xold, double xnew) const
{
    const double dold = std::floor(xold), dnew = std::floor(xnew);
    return std::log((raw_cdf(dnew + 1.) - raw_cdf(dnew)) / (raw_cdf(dold + 1.) - raw_cdf(dold)));
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

double Distribution::__trunc_dummy(double x) const { return x; }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  Base Distribution : truncation (continuos)  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__trunc_min(double x) const
{
    return x / (1. - raw_cdf(get_min()));
}
double Distribution::__trunc_max(double x) const
{
    return x / (raw_cdf(get_max()));
}
double Distribution::__trunc_both_sides(double x) const
{
    return x / (raw_cdf(get_max()) - raw_cdf(get_min()));
}
double Distribution::__trunc_rand(rnd_gen &prng)
{
    return _rejection_sampling(prng, &Distribution::__default_rand);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////  Base Distribution : truncation (discrete defined)  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__trunc_min_discr_defined(double x) const
{
    return x / (1. - raw_cdf(get_min() - 1));
}
double Distribution::__trunc_max_discr_defined(double x) const
{
    return x / (raw_cdf(get_max() + 1));
}
double Distribution::__trunc_both_sides_discr_defined(double x) const
{
    return x / (raw_cdf(get_max() + 1) - raw_cdf(get_min() - 1));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  Base Distribution : truncation (CENTER discretization)  ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//----- Case of CENTER Discretization -----//
double Distribution::__trunc_min_discr_center(double x) const
{
    return x / (1. - raw_cdf(std::ceil(get_min()) - .5));
}
double Distribution::__trunc_max_discr_center(double x) const
{
    return x / (raw_cdf(std::floor(get_max()) + .5));
}
double Distribution::__trunc_both_sides_discr_center(double x) const
{
    return x / (raw_cdf(std::floor(get_max()) + .5) - raw_cdf(std::ceil(get_min()) - .5));
}
double Distribution::__trunc_rand_discr_center(rnd_gen &prng)
{
    return _rejection_sampling(prng, &Distribution::__discr_center_rand);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Base Distribution : truncation (LEFT discretization)  /////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__trunc_min_discr_left(double x) const
{
    return x / (1. - raw_cdf(std::floor(get_min())));
}
double Distribution::__trunc_max_discr_left(double x) const
{
    return x / (raw_cdf(std::ceil(get_max() + 1)));
}
double Distribution::__trunc_both_sides_discr_left(double x) const
{
    return x / (raw_cdf(std::ceil(get_max() + 1)) - raw_cdf(std::floor(get_min())));
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

double Distribution::__logtrunc_dummy(double x) const
{
    return x;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Base Distribution : log-truncation (continuos)  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_min(double x) const
{
    return x - std::log(1. - raw_cdf(get_min()));
}
double Distribution::__logtrunc_max(double x) const
{
    return x - std::log(raw_cdf(get_max()));
}
double Distribution::__logtrunc_both_sides(double x) const
{
    return x - std::log(raw_cdf(get_max()) - raw_cdf(get_min()));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////  Base Distribution : log-truncation (discrete defined)  ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_min_discr_defined(double x) const
{
    return x - std::log(1. - raw_cdf(get_min() - 1));
}
double Distribution::__logtrunc_max_discr_defined(double x) const
{
    return x - std::log(raw_cdf(get_max() + 1));
}
double Distribution::__logtrunc_both_sides_discr_defined(double x) const
{
    return x - std::log(raw_cdf(get_max() + 1) - raw_cdf(get_min() - 1));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////  Base Distribution : log-truncation (CENTER discretization)  /////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_min_discr_center(double x) const
{
    return x - std::log(1. - raw_cdf(std::ceil(get_min()) - .5));
}
double Distribution::__logtrunc_max_discr_center(double x) const
{
    return x - std::log(raw_cdf(std::floor(get_max()) + .5));
}
double Distribution::__logtrunc_both_sides_discr_center(double x) const
{
    return x - std::log(raw_cdf(std::floor(get_max()) + .5) - raw_cdf(std::ceil(get_min()) - .5));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  Base Distribution : log-truncation (LEFT discretization)  //////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::__logtrunc_min_discr_left(double x) const
{
    return x - std::log(1. - raw_cdf(std::floor(get_min())));
}
double Distribution::__logtrunc_max_discr_left(double x) const
{
    return x - std::log(raw_cdf(std::ceil(get_max())));
}
double Distribution::__logtrunc_both_sides_discr_left(double x) const
{
    return x - std::log(raw_cdf(std::ceil(get_max())) - raw_cdf(std::floor(get_min())));
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

double Distribution::__shift_dummy(double x) const
{
    return x;
}
double Distribution::__shift_min(double x) const
{
    return x - raw_cdf(get_min());
}
double Distribution::__shift_min_discr_defined(double x) const
{
    return x - raw_cdf(get_min() - 1);
}
double Distribution::__shift_min_discr_center(double x) const
{
    return x - raw_cdf(std::floor(get_min()) - .5);
}
double Distribution::__shift_min_discr_left(double x) const
{
    return x - raw_cdf(std::floor(get_min()));
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
bool Distribution::is_discrete() const { return _do_pdf != &Distribution::__default_pdf; }

/* Select functional method in case of discretization. */
void Distribution::init_functional_discretization(ushort discr_type)
{
    //===== DISCRETIZATION =====//
    switch (discr_type)
    {
    case D_DISCR_TYPE_DFND: // discrete defined
        _do_pdf = &Distribution::__discr_defined_pdf;
        _do_cdf = &Distribution::__discr_defined_cdf;
        _do_logpdf = &Distribution::__discr_defined_logpdf;
        _do_logratio = &Distribution::__discr_defined_logratio;
        _do_rand = &Distribution::__default_rand;
        break;

    case D_DISCR_TYPE_CNTR: // centered discretization
        _do_pdf = &Distribution::__discr_center_pdf;
        _do_cdf = &Distribution::__discr_center_cdf;
        _do_logpdf = &Distribution::__discr_center_logpdf;
        _do_logratio = &Distribution::__discr_center_logratio;
        _do_rand = &Distribution::__discr_center_rand;
        break;

    case D_DISCR_TYPE_LEFT: // left discretization
        _do_pdf = &Distribution::__discr_left_pdf;
        _do_cdf = &Distribution::__discr_left_cdf;
        _do_logpdf = &Distribution::__discr_left_logpdf;
        _do_logratio = &Distribution::__discr_left_logratio;
        _do_rand = &Distribution::__discr_left_rand;
        break;

    default: // no discretization (default)
        _do_pdf = &Distribution::__default_pdf;
        _do_cdf = &Distribution::__default_cdf;
        _do_logpdf = &Distribution::__default_logpdf;
        _do_logratio = &Distribution::__default_logratio;
        _do_rand = &Distribution::__default_rand;
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
            _do_trunc = &Distribution::__trunc_min_discr_defined;
            _do_logtrunc = &Distribution::__logtrunc_min_discr_defined;
            _do_shift = &Distribution::__shift_min_discr_defined;
            _do_rand = &Distribution::__trunc_rand;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            _do_trunc = &Distribution::__trunc_min_discr_center;
            _do_logtrunc = &Distribution::__logtrunc_min_discr_center;
            _do_shift = &Distribution::__shift_min_discr_center;
            _do_rand = &Distribution::__trunc_rand_discr_center;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _do_trunc = &Distribution::__trunc_min_discr_left;
            _do_logtrunc = &Distribution::__logtrunc_min_discr_left;
            _do_shift = &Distribution::__shift_min_discr_left;
            _do_rand = &Distribution::__trunc_rand_discr_left;
            break;
        default: // no discretization
            _do_trunc = &Distribution::__trunc_min;
            _do_logtrunc = &Distribution::__logtrunc_min;
            _do_shift = &Distribution::__shift_min;
            _do_rand = &Distribution::__trunc_rand;
        }
        break;

    case 2: //----- upper range truncation -----//
        _do_shift = &Distribution::__shift_dummy;
        switch (discr_type)
        {
        case D_DISCR_TYPE_DFND: // discrete defined
            _do_trunc = &Distribution::__trunc_max_discr_defined;
            _do_logtrunc = &Distribution::__logtrunc_max_discr_defined;
            _do_rand = &Distribution::__trunc_rand;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            _do_trunc = &Distribution::__trunc_max_discr_center;
            _do_logtrunc = &Distribution::__logtrunc_max_discr_center;
            _do_rand = &Distribution::__trunc_rand_discr_center;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _do_trunc = &Distribution::__trunc_max_discr_left;
            _do_logtrunc = &Distribution::__logtrunc_max_discr_left;
            _do_rand = &Distribution::__trunc_rand_discr_left;
            break;
        default: // no discretization
            _do_trunc = &Distribution::__trunc_max;
            _do_logtrunc = &Distribution::__logtrunc_max;
            _do_rand = &Distribution::__trunc_rand;
        }
        break;

    case 3: //----- both sides truncation -----//
        switch (discr_type)
        {
        case D_DISCR_TYPE_DFND: // discrete defined
            _do_trunc = &Distribution::__trunc_both_sides_discr_defined;
            _do_logtrunc = &Distribution::__logtrunc_both_sides_discr_defined;
            _do_shift = &Distribution::__shift_min_discr_defined;
            _do_rand = &Distribution::__trunc_rand;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            _do_trunc = &Distribution::__trunc_both_sides_discr_center;
            _do_logtrunc = &Distribution::__logtrunc_both_sides_discr_center;
            _do_shift = &Distribution::__shift_min_discr_center;
            _do_rand = &Distribution::__trunc_rand_discr_center;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            _do_trunc = &Distribution::__trunc_both_sides_discr_left;
            _do_logtrunc = &Distribution::__logtrunc_both_sides_discr_left;
            _do_shift = &Distribution::__shift_min_discr_left;
            _do_rand = &Distribution::__trunc_rand_discr_left;
            break;
        default: // no discretization
            _do_trunc = &Distribution::__trunc_both_sides;
            _do_logtrunc = &Distribution::__logtrunc_both_sides;
            _do_shift = &Distribution::__shift_min;
            _do_rand = &Distribution::__trunc_rand;
        }
        break;

    default: //----- no truncation (case 0) -----//
        _do_trunc = &Distribution::__trunc_dummy;
        _do_logtrunc = &Distribution::__logtrunc_dummy;
        _do_shift = &Distribution::__shift_dummy;
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
const std::vector<dbl_vec> Distribution::showcase(int n_bins, double ytol) const
{
    if (n_bins < 1)
        throw std::invalid_argument("Number of bins should be greater than 0.");

    // pick extremes for the showcase
    double lower = get_min();
    double upper = get_max();

    if (tools::is_na(upper))
    {
        upper = 10. + (tools::is_na(lower) ? 0 : lower);
        while ((1. - cdf(upper) > ytol) && (upper < DEFAULT_MAX_SIZE))
            ++upper;
    }

    if (tools::is_na(lower))
    {
        lower = -10.;
        while ((cdf(lower) > ytol) && (lower > -DEFAULT_MAX_SIZE))
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
        inst[i] = pdf(x);
        cumul[i] = cdf(x);
    }

    return {gentime, inst, cumul};
}

/* Developer method to check for normalization. */
double Distribution::__d_e_v__check__pdf__norm(int n_bins,
                                               double ytol) const
{
    // produce a showcase of the distribution profile
    const std::vector<dbl_vec> values = showcase(n_bins, ytol);

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
