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
        return (this->*logratio_)(xold, xnew);
}
double Distribution::rand(rnd_gen &prng)
{
    update_args(); // UPDATE: random generator classes should be redisgned to deal with this
    return (this->*rand_)(prng);
}

/* Function-like return value for probability density. */
double Distribution::do_return_pdf_(double x) const
{
    if (valid(x))
        return pdf(x);
    else
        return 0.;
}

/* Function-like return value for cumulative distribution. */
double Distribution::do_return_cdf_(double x) const
{
    if (below(x))
        return 0.;
    else if (above(x))
        return cdf(get_max());
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
double Distribution::rejection_sampling_(rnd_gen &prng,
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

double Distribution::do_discr_center_pdf_(double x) const
{
    const double d = round(x);
    return (this->*trunc_)(raw_cdf(d + .5) - raw_cdf(d - .5));
}

double Distribution::do_discr_center_cdf_(double x) const
{
    const double d = round(x);
    return (this->*trunc_)((this->*shift_)(raw_cdf(d + .5)));
}

double Distribution::do_discr_center_logpdf_(double x) const
{
    return std::log(do_discr_center_pdf_(x));
}

double Distribution::do_discr_center_logratio_(double xold, double xnew) const
{
    const double dold = std::round(xold), dnew = std::round(xnew);
    return std::log((raw_cdf(dnew + .5) - raw_cdf(dnew - .5)) / (raw_cdf(dold + .5) - raw_cdf(dold - .5)));
}

double Distribution::do_discr_center_rand_(rnd_gen &prng) { return static_cast<double>(std::lround(do_default_rand_(prng))); }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Base Distribution : discretized LEFT  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

double Distribution::do_discr_left_pdf_(double x) const
{
    const double d = std::floor(x);
    return (this->*trunc_)(raw_cdf(d + 1.) - raw_cdf(d));
}

double Distribution::do_discr_left_cdf_(double x) const
{
    return (this->*trunc_)((this->*shift_)(raw_cdf(std::floor(x) + 1.)));
}

double Distribution::do_discr_left_logpdf_(double x) const
{
    return std::log(do_discr_left_pdf_(x));
}

double Distribution::do_discr_left_logratio_(double xold, double xnew) const
{
    const double dold = std::floor(xold), dnew = std::floor(xnew);
    return std::log((raw_cdf(dnew + 1.) - raw_cdf(dnew)) / (raw_cdf(dold + 1.) - raw_cdf(dold)));
}

double Distribution::do_discr_left_rand_(rnd_gen &prng)
{
    return std::floor(do_default_rand_(prng));
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
        pdf_ = &Distribution::do_discr_defined_pdf_;
        cdf_ = &Distribution::do_discr_defined_cdf_;
        logpdf_ = &Distribution::do_discr_defined_logpdf_;
        logratio_ = &Distribution::do_discr_defined_logratio_;
        rand_ = &Distribution::do_default_rand_;
        break;

    case D_DISCR_TYPE_CNTR: // centered discretization
        pdf_ = &Distribution::do_discr_center_pdf_;
        cdf_ = &Distribution::do_discr_center_cdf_;
        logpdf_ = &Distribution::do_discr_center_logpdf_;
        logratio_ = &Distribution::do_discr_center_logratio_;
        rand_ = &Distribution::do_discr_center_rand_;
        break;

    case D_DISCR_TYPE_LEFT: // left discretization
        pdf_ = &Distribution::do_discr_left_pdf_;
        cdf_ = &Distribution::do_discr_left_cdf_;
        logpdf_ = &Distribution::do_discr_left_logpdf_;
        logratio_ = &Distribution::do_discr_left_logratio_;
        rand_ = &Distribution::do_discr_left_rand_;
        break;

    default: // no discretization (default)
        pdf_ = &Distribution::do_default_pdf_;
        cdf_ = &Distribution::do_default_cdf_;
        logpdf_ = &Distribution::do_default_logpdf_;
        logratio_ = &Distribution::do_default_logratio_;
        rand_ = &Distribution::do_default_rand_;
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
            trunc_ = &Distribution::do__trunc_min_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            trunc_ = &Distribution::do_trunc_min_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_min_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_left_;
            shift_ = &Distribution::do_shift_min_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_min_;
            logtrunc_ = &Distribution::do_logtrunc_min_;
            shift_ = &Distribution::do_shift_min_;
            rand_ = &Distribution::do_trunc_rand_;
        }
        break;

    case 2: //----- upper range truncation -----//
        shift_ = &Distribution::do_shift_dummy_;
        switch (discr_type)
        {
        case D_DISCR_TYPE_DFND: // discrete defined
            trunc_ = &Distribution::do_trunc_max_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            trunc_ = &Distribution::do_trunc_max_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_max_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_max_;
            logtrunc_ = &Distribution::do_logtrunc_max_;
            rand_ = &Distribution::do_trunc_rand_;
        }
        break;

    case 3: //----- both sides truncation -----//
        switch (discr_type)
        {
        case D_DISCR_TYPE_DFND: // discrete defined
            trunc_ = &Distribution::do_trunc_both_sides_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case D_DISCR_TYPE_CNTR: // centered discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case D_DISCR_TYPE_LEFT: // left discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_left_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_left_;
            shift_ = &Distribution::do_shift_min_discr_left_;
            rand_ = &Distribution::do_trunc_rand_discr_left_;
            break;
        default: // no discretization
            trunc_ = &Distribution::do_trunc_both_sides_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_;
            shift_ = &Distribution::do_shift_min_;
            rand_ = &Distribution::do_trunc_rand_;
        }
        break;

    default: //----- no truncation (case 0) -----//
        trunc_ = &Distribution::do_trunc_dummy_;
        logtrunc_ = &Distribution::do_logtrunc_dummy_;
        shift_ = &Distribution::do_shift_dummy_;
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
double Distribution::develop_check_pdf_norm_(int n_bins, double ytol) const
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

const str_vec DistributionRegistry::get_options() const
{
    str_vec keys;
    keys.reserve(_factory_map.size());
    for (const auto &pair : _factory_map)
        keys.push_back(pair.first);
    return keys;
}

/* Try for default distribution name, otherwise with R-like prefix. */
bool DistributionRegistry::is_valid_option(const std::string &opt) const
{
    const str_vec options = get_options();
    if (options % opt)
        return true;
    else if (options % opt.substr(1))
        return str_vec({"d", "p"}) % opt.substr(0, 1);
    return false;
}

std::shared_ptr<Function> DistributionRegistry::get(const std::string &opt) const
{
    const str_vec options = get_options();
    bool is_R_like_prefix = false;

    // try for default name or R-like prefix
    if (!(options % opt))
        if (options % opt.substr(1))
            is_R_like_prefix = true;

    auto output = is_R_like_prefix ? get_distr(opt.substr(1)) : get_distr(opt);

    // if there's a prefix, see if cumulative and in case set it
    if (is_R_like_prefix)
        if (opt.substr(0, 1) == "p")
            output->set_return_cumul();

    return output;
}

std::shared_ptr<Distribution> DistributionRegistry::get_distr(const std::string &opt) const
{
    if (_factory_map.find(opt) != _factory_map.end())
        return _factory_map.find(opt)->second();

    else if (_factory_map.find(opt.substr(1)) != _factory_map.end())
        return _factory_map.find(opt.substr(1))->second();

    return nullptr;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Distribution Registry : setter  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<Distribution> compose_distr_registries(const std::string &opt,
                                                       const std::string &name,
                                                       const std::vector<const DistributionRegistry *> &registries)
{
    str_vec options = {};

    if (registries.empty())
        throw std::invalid_argument("No registries provided for \"" + name + "\".");
    else if (tools::is_na(opt))
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
