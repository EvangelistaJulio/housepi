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

double Distribution::logratio(double xold, double xnew) const
{
    if (xold == xnew)
        return 0.;
    else
        return (this->*_do_logratio)(xold, xnew);
}
double Distribution::rand(rnd_gen &prng)
{
    update_args(); // UPDATE: random generator classes should be redisgned to deal with this
    return (this->*_do_rand)(prng);
}

/* Function-like return value for probability density. */
double Distribution::__return_pdf(double x) const
{
    if (valid(x))
        return pdf(x);
    else
        return 0.;
}

/* Function-like return value for cumulative distribution. */
double Distribution::__return_cdf(double x) const
{
    if (below(x))
        return 0.;
    else if (above(x))
        return 1.;
    else
        return cdf(x);
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
///////////////////////////  Base Distribution : functional reshaping  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

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
    init_functional_range(false, false);
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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Distribution Registry : getter  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

const str_vec DistrRegistry::get_options() const
{
    str_vec keys;
    keys.reserve(_factory_map.size());
    for (const auto &pair : _factory_map)
        keys.push_back(pair.first);
    return keys;
}

/* Try for default distribution name, otherwise with R-like prefix. */
bool DistrRegistry::is_valid_option(const std::string &opt) const
{
    const str_vec options = get_options();
    if (contains(options, opt))
        return true;
    else if (contains(options, opt.substr(1)))
        return true;
    return false;
}

std::shared_ptr<Distribution> DistrRegistry::get_distr(const std::string &opt) const
{
    auto it = _factory_map.find(opt);
    if (it != _factory_map.end())
        return it->second();
    return nullptr;
}

std::shared_ptr<Function> DistrRegistry::get(const std::string &opt) const
{
    const str_vec options = get_options();
    bool is_R_like_prefix = false;

    // try for default name or R-like prefix
    if (!contains(options, opt))
        if (contains(options, opt.substr(1)))
            is_R_like_prefix = true;

    auto output = is_R_like_prefix ? get_distr(opt.substr(1)) : get_distr(opt);

    // if there's a prefix, see if cumulative and in case set it
    if (is_R_like_prefix)
        if (opt.substr(0, 1) == "p")
            output->set_return_cumul();

    return output;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Distribution Registry : setter  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void DistrRegistry::register_distr(const std::string &opt, FactoryDistr creator)
{
    _factory_map[opt] = creator;
}

void DistrRegistry::register_func(const std::string &, FactoryFunc) { throw std::invalid_argument(
    "Impossible to register a function in distribution registry \"" + name() + "\"."); }

std::shared_ptr<Distribution> compose_distr_registries(const std::string &opt,
                                                       const std::string &name,
                                                       const std::vector<const DistrRegistry *> &registries)
{
    str_vec options = {};

    if (tools::is_na(opt))
        return std::make_shared<Distribution>();
    else
        for (auto &r : registries)
            if (r->is_valid_option(opt))
                return r->get_distr(opt);
            else
                options += r->get_options();

    throw std::invalid_argument(
        "Unknown " + name + " \"" + opt + "\"." + tools::find_closest(opt, options));
}
