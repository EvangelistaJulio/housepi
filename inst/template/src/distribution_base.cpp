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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Distribution - rand  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// UPDATE: random generator classes should be redesigned to deal with this
double Distribution::rand(rnd_gen &prng)
{
    update_args_();
    return (this->*rand_)(prng);
}

/// Rejection sampling.
double Distribution::rejec_sample_(rnd_gen &prng, double (Distribution::*func)(rnd_gen &))
{
    double xnew = (this->*func)(prng);

    unsigned n_attempts = 0;
    while (!is_valid(xnew))
    {
        if (n_attempts > HSP_MAX_NUM_ATTEMPTS)
            throw std::runtime_error(
                "Exceeded attempts (N > " + std::to_string(HSP_MAX_NUM_ATTEMPTS) + ") in random draw of " + str_summary() + ".");
        xnew = (this->*func)(prng);
        ++n_attempts;
    }

    return xnew;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Distribution - reshaping  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Select functional for discretization of random sampling.
/// UPDATE: to be absorbed inside init_functional_discretization_
void Distribution::init_sampling_discretization_()
{
    unsigned short discr_type = which_discrete_case_();

    switch (discr_type)
    {
    case HSP_DISCR_DEFINED: // discrete defined
        rand_ = &Distribution::do_rand_;
        break;

    case HSP_DISCR_CENTER: // centered discretization
        rand_ = &Distribution::do_discr_center_rand_;
        break;

    case HSP_DISCR_LEFT: // left discretization
        rand_ = &Distribution::do_discr_left_rand_;
        break;

    default: // no discretization (default)
        rand_ = &Distribution::do_rand_;
    }
}

/// Select functional method in case of truncation.
void Distribution::init_functional_truncation_()
{
    unsigned short discr_type = which_discrete_case_();

    switch (which_valid_case_())
    {
    case 1: //----- lower range truncation -----//
        switch (discr_type)
        {
        case HSP_DISCR_DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_min_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case HSP_DISCR_CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_min_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_min_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case HSP_DISCR_LEFT: // left discretization
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
        case HSP_DISCR_DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_max_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case HSP_DISCR_CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_max_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_max_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case HSP_DISCR_LEFT: // left discretization
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
        case HSP_DISCR_DEFINED: // discrete defined
            trunc_ = &Distribution::do_trunc_both_sides_discr_defined_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_defined_;
            shift_ = &Distribution::do_shift_min_discr_defined_;
            rand_ = &Distribution::do_trunc_rand_;
            break;
        case HSP_DISCR_CENTER: // centered discretization
            trunc_ = &Distribution::do_trunc_both_sides_discr_center_;
            logtrunc_ = &Distribution::do_logtrunc_both_sides_discr_center_;
            shift_ = &Distribution::do_shift_min_discr_center_;
            rand_ = &Distribution::do_trunc_rand_discr_center_;
            break;
        case HSP_DISCR_LEFT: // left discretization
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

/// Default reshaping initialization.
/// note: Necessary to have dynamic state variables for NULL case
void Distribution::init_reshaping(bool discrete)
{
    init_functional_discretization_((discrete) ? HSP_DISCR_DEFINED : HSP_DISCR_NONE);
    init_sampling_discretization_();
    init_functional_range_(false, false);
    init_functional_truncation_();
}

/// Functional range for lower bound and no upper bound.
bool Distribution::init_range_min_(double min_val, bool is_open_zero)
{
    if (get_min() < min_val)
        throw std::invalid_argument(
            "Invalid min value " + str_min() + " for distribution \"" + name() + "\".");
    else if ((get_min() > min_val) || (!tools::is_na(get_max())))
    {
        init_functional_range_(((is_open_zero) ? (get_min() == min_val) : false), false);
        return true;
    }
    return false;
}

/// Frequent reshaping: domain has zero lower bound.
void Distribution::init_reshaping_min_zero_(bool discrete, bool is_open_zero)
{
    init_functional_discretization_((discrete) ? HSP_DISCR_LEFT : HSP_DISCR_NONE);
    init_sampling_discretization_();
    if (init_range_min_(0., is_open_zero))
        init_functional_truncation_();
}

/// Frequent reshaping: domain has zero lower bound for discrete defined distributions.
void Distribution::init_reshaping_min_zero_discr_defined_()
{
    init_functional_discretization_(HSP_DISCR_DEFINED);
    init_sampling_discretization_();
    if (init_range_min_(0., false))
        init_functional_truncation_();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Distribution Registry : getter  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Override get_options method of generic function with distribution.
/// note: this differ in the choice of factory_map_ member.
str_vec DistributionRegistry::get_options() const
{
    str_vec keys;
    keys.reserve(factory_map_.size());
    for (const auto &pair : factory_map_)
        keys.push_back(pair.first);
    return keys;
}

std::shared_ptr<Distribution> DistributionRegistry::get_distr(const std::string &opt) const
{
    if (factory_map_.find(opt) != factory_map_.end())
        return factory_map_.find(opt)->second();

    return nullptr;
}

std::shared_ptr<Function> DistributionRegistry::get(const std::string &opt, bool is_primitive) const
{
    auto output = get_distr(opt);
    if (is_primitive)
        output->set_return_as_distr_primitive();
    else
        output->set_return_as_distr_call();
    return output;
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
        tools::str::find_closest(opt, options, HSP_MIN_STR_DIST));
}
