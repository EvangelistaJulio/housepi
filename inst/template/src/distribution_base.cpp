/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
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

double Distribution::rand(rnd_gen &prng)
{
    update_args(); // UPDATE: random generator classes should be redesigned to deal with this
    return (this->*rand_)(prng);
}
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Rejection sampling. */
double Distribution::rejec_sample_(rnd_gen &prng, double (Distribution::*func)(rnd_gen &))
{
    double xnew = (this->*func)(prng);

    unsigned n_attempts = 0;
    while (!is_valid(xnew))
    {
        if (n_attempts > MAX_NUM_ATTEMPTS)
            throw std::runtime_error(
                "Exceeded attempts (N > " + std::to_string(MAX_NUM_ATTEMPTS) + ") in random draw of " + str_summary() + ".");
        xnew = (this->*func)(prng);
        ++n_attempts;
    }

    return xnew;
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
void Distribution::init_functional_discretization(unsigned short discr_type)
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
void Distribution::init_functional_truncation(unsigned short discr_type)
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
    unsigned short discr_case = (discrete) ? D_DISCR_TYPE_DFND : D_DISCR_TYPE_NONE;
    init_functional_discretization(discr_case);
    init_functional_range(false, false);
    init_functional_truncation(discr_case);
}

/* Functional range for lower bound and no upper bound. */
bool Distribution::init_range_min(double min_val, bool is_open_zero)
{
    if (get_min() < min_val)
        throw std::invalid_argument(
            "Invalid min value " + str_min() + " for distribution \"" + name() + "\".");
    else if ((get_min() > min_val) || (!tools::is_na(get_max())))
    {
        init_functional_range(((is_open_zero) ? (get_min() == min_val) : false), false);
        return true;
    }
    return false;
}

/* Frequent reshaping: domain has zero lower bound. */
void Distribution::init_reshaping_min_zero(bool discrete, bool is_open_zero)
{
    unsigned short discr_case = (discrete) ? D_DISCR_TYPE_LEFT : D_DISCR_TYPE_NONE;
    init_functional_discretization(discr_case);
    if (init_range_min(0., is_open_zero))
        init_functional_truncation(discr_case);
}

/* Frequent reshaping: domain has zero lower bound for discrete defined distributions. */
void Distribution::init_reshaping_min_zero_discr_defined()
{
    init_functional_discretization(D_DISCR_TYPE_DFND);
    if (init_range_min(0., false))
        init_functional_truncation(D_DISCR_TYPE_DFND);
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

str_vec DistributionRegistry::get_options() const
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
        "Unknown " + name + " \"" + opt + "\"." + 
        tools::str::find_closest(opt, options, MIN_STR_DIST));
}
