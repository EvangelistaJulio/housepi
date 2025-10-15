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
 * mcmc_build_param.cpp
 *
 */

#include <atomic>

#include "tools/error.h"
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC::Impl : evaluate loglikelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize model and compute the loglikelihood.
double MCMC::Impl::evaluate()
{
    try
    {
        setup.n_chains = 0;          // use 0 chains as a state variable
        assert_proper_dataset_();    // check if dataset is empty
        param.user_reset();          // initialize non constant parameters
        augm.init();                 // initialize individual augmented relations
        augm.user_reset();           // initialize non constant augmented
        init_model_for_inference_(); // initialize model for inference
        update_full_loglike_(true);  // initialize cache + full loglikelihood
        return full_loglike_;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error(e.what());
    }

    return hsp::na::DOUBLE;
}

/// Assert dataset is not empty.
void MCMC::Impl::assert_proper_dataset_()
{
    if (dataset.n_houses() == 0)
        throw std::length_error("Unable to proceed: no data was found.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : update house loglikelihood  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute new loglikelihood, store old value in buffer. Same for households.
void MCMC::Impl::update_full_loglike_(bool cache_depen)
{
    const bool is_augm_ll = augm.is_any_model();
    double new_ll = 0.0;
    std::atomic<bool> has_error(false);

    // buffer infectivity profile if lookup tables are active
    infec_profile.update(cache_depen);

#pragma omp parallel for reduction(+ : new_ll)
    for (size_t R = 0; R < dataset.n_restricted_houses(); ++R)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        Household &this_house = dataset.restricted(R);
        try
        {
            new_ll += update_house_loglike_(this_house, is_augm_ll, cache_depen);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    error_message = "While computing loglikelihood for " +
                                    this_house.err_name() + ":\n" +
                                    std::string(e.what());
            }
        }
    }
    if (!error_message.empty())
        throw std::runtime_error(error_message);

    set_full_loglike_(new_ll);
}

/// Update and return contribution of given household to full logliklehood.
/// NOTE: house multiplicity doesn't work with data augmentation
double MCMC::Impl::update_house_loglike_(Household &house, bool is_augm_ll, bool is_cache_depen)
{
    // buffer all household cache dependencies
    if (is_cache_depen)
    {
        house.buffer_all_cache();
        for (Individual &member : house.infected())
            update_cache_dependecies_(house, member);
    }

    // update household observation loglikelihood
    double output = compute_and_set_household_obs_loglike_(house);

    // update individual augmented loglikelihood
    if (is_augm_ll)
        for (Individual &member : house.infected())
            output += member.compute_and_set_augm_model_loglike();

    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : Impl : inference  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Performs the inference using Metropolis-Hastings.
void MCMC::Impl::infer()
{
    try
    {
        setup.prog_declare();                               // declare current progress bar
        init_monte_carlo_();                                // initialize mcmc sampling
        init_model_for_inference_();                        // user initialization model
        augm.init();                                        // initialize augmented variables
        setup.prog_open(setup.c_n_sweeps * setup.n_chains); // open progress bar
        run_monte_carlo_();                                 // run Monte Carlo
        setup.prog_close();                                 // close progress bar
    }
    catch (std::exception &e)
    {
        setup.prog_crash();
        hsp::msg::error(e.what());
    }
}

/// Assert setup definitions are ok to proceed with mcmc sampling.
void MCMC::Impl::init_monte_carlo_()
{
    setup.assert_mcmc_setup();

    // check if there is a dataset
    assert_proper_dataset_();

    // check if there is actually something to sample
    if (!(param.is_any_sample() || augm.is_any_sample()))
        throw std::length_error(
            "Nothing to sample. Check your parameter and augmented definitions.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : run monte carlo  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Markov chain core, used for inference and power posteriors.
/// NOTE: sweep is the number of previously completed thin sweeps.
void MCMC::Impl::run_monte_carlo_()
{
    unsigned long &chain = setup.chain; // current chain
    init_posterior_sample_();           // initialize diagnostics

    // run multiple MCMC chains&
    for (chain = 1; chain <= setup.n_chains; ++chain)
    {
        setup.prog_set_text(); // update progress bar text
        param.user_reset();    // reset parameter user definitions
        augm.user_reset();     // reset augmented user definitions
        update_full_loglike_(true); // init cache + full loglikelihood
        init_move_counts_();        // init move counts for adaptive

        // (Finite) Adaptative Iterations
        for (unsigned long i = 0; i < setup.n_adapt; i += setup.n_thin, setup.prog_update())
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_and_adapt_params_(i + t);

        // Burn-in Iterations
        for (unsigned long i = 0; i < setup.n_burnin; i += setup.n_thin, setup.prog_update())
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_();

        // Posterior sampling
        for (unsigned long i = setup.n_burnin; i < setup.n_steps; i += setup.n_thin, setup.prog_update())
        {
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_();
            update_posterior_sample_(); // store sample
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : monte carlo moves ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Move all model parameters.
void MCMC::Impl::move_all_()
{
    // move augmented variables if any
    move_all_augm_();

    // move all simple model parameters
    for (const std::string &p : param.list_sample_standard)
        move_single_parameter_(p, false);

    // move parameters with cache dependecies
    for (const std::string &p : param.list_sample_cache)
        move_single_parameter_(p, true);
}

/// Attempt to move a single parameter for the given inverse temperature.
bool MCMC::Impl::move_single_parameter_(const std::string &p, bool cache_depen)
{
    // proceed only if the new value is consistent
    // otherwise automatically refuse move
    if (propose_valid_param_move_(p))
    {
        // if parameter is cache dependent, update cache
        update_full_loglike_(cache_depen);

        // perform Metropolis-Hastings to accept or not move
        if (accept_valid_param_move_(p))
            return true; // accepted move
        else
            recover_parameter_(p, cache_depen);
    }
    return false; // refused move
}

/// Move all augmented variables and model parameters and adapt their rule arguments.
void MCMC::Impl::move_all_and_adapt_params_(unsigned long step)
{
    move_all_augm_(); // move augmented variables to sample

    // move all model parameters and adapt their rule arguments
    const double ratio = static_cast<double>(step) / setup.n_adapt;

    // UPDATE: loops can be improved with viewers

    // move all simple model parameters
    for (const std::string &p : param.list_sample_standard)
    {
        Variable &this_param = param.at(p);
        this_param.update_accept_rate(move_single_parameter_(p, false));
        this_param.adapt_rule(ratio);
    }

    // move parameters with cache dependecies
    for (const std::string &p : param.list_sample_cache)
    {
        Variable &this_param = param.at(p);
        this_param.update_accept_rate(move_single_parameter_(p, true));
        this_param.adapt_rule(ratio);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : parameter moves  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Set accept rates to zero.
void MCMC::Impl::init_move_counts_()
{
    for (const std::string &p : param.list_sample)
        param.at(p).reset_accept_rate();
}

/// propose a valid move (return true), otherwise update refusal (return false).
bool MCMC::Impl::propose_valid_param_move_(const std::string &p)
{
    Variable &this_par = param.at(p);

    // set new value of the parameter and store the old to buffer
    this_par.propose(setup.prng());

    // recover buffer if not valid move
    if (!this_par.is_valid_move())
    {
        this_par.recover_value();
        return false;
    }

    return true;
}

/// Restore the parameter value in the buffer.
void MCMC::Impl::recover_parameter_(const std::string &p, bool is_cache_depen)
{
    // recover buffer value of parameter
    param.at(p).recover_value();

    // recover buffer value of full loglike
    full_loglike_ = full_buffer_;

    // recover buffered loglikelihood in households
    for (Household &house : dataset.restricted())
    {
        house.recover_buffer_values(is_cache_depen);
        if (augm.is_any_sample())
            for (Individual &member : house.infected())
                member.recover_augm_model_loglike();
    }

    infec_profile.recover(is_cache_depen);
}

/// Metropolis-Hastings to accept or not current parameter move.
bool MCMC::Impl::accept_valid_param_move_(const std::string &p)
{
    const double D_ll = full_loglike_ - full_buffer_;
    const double move_logratio = param.at(p).get_move_logratio();
    return setup.unif() < std::exp(setup.inv_T * D_ll + move_logratio);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : posterior sample  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize posterior sample of parameters.
void MCMC::Impl::Parameter::init_draws(unsigned long sample_size)
{
    smpl_values.clear();
    smpl_values.reserve(list_sample.size());
    for (const std::string &p : list_sample)
    {
        smpl_values[p].clear();
        smpl_values[p].reserve(sample_size);
        smpl_values[p].assign(sample_size, hsp::na::DOUBLE);
    }
}

/// Store current sample of parameters.
void MCMC::Impl::Parameter::update_draws(unsigned long smpl_count)
{
    for (const std::string &p : list_sample)
        smpl_values[p][smpl_count] = get(p);
}

/// Deep clear existing posterior sample of parameters.
void MCMC::Impl::Parameter::clear_draws()
{
    for (auto &pair : smpl_values)
    {
        pair.second.clear();
        pair.second.shrink_to_fit();
    }
    smpl_values.clear();
    smpl_values.rehash(0);
}

/// Initialize posterior sample of parameters.
void MCMC::Impl::DataWide::init_draws(unsigned long sample_size)
{
    // sample over households
    smpl_values_.clear();
    smpl_values_.resize(restrict_indexes_.size());
    for (size_t R = 0; R < restrict_indexes_.size(); ++R)
    {
        smpl_values_[R].clear();
        smpl_values_[R].reserve(sample_size);
        smpl_values_[R].assign(sample_size, hsp::na::DOUBLE);
    }
}

/// Store current sample of loglikelihood for restricted households.
void MCMC::Impl::DataWide::update_draws(unsigned long smpl_count)
{
    // store sample for restricted households
    for (size_t R = 0; R < restrict_indexes_.size(); ++R)
    {
        Household &house = house_vector_[restrict_indexes_[R]];
        double x = house.get_obs_loglike() / house.multiplicity();
        for (const Individual &member : house.infected())
            x += member.get_augm_model_loglike();
        smpl_values_[R][smpl_count] = x;
    }
}

/// Deep clear posterior sample of pointwise loglikelihood.
void MCMC::Impl::DataWide::clear_draws()
{
    for (auto &vec : smpl_values_)
    {
        vec.clear();
        vec.shrink_to_fit();
    }
    smpl_values_.clear();
    smpl_values_.shrink_to_fit();
}

/// Init posterior sampling.
void MCMC::Impl::init_posterior_sample_()
{
    unsigned long sample_size = setup.c_sample_size * setup.n_chains;

    // diagnostics over full loglikelihood
    smpl_count_ = 0;
    smpl_full_loglike_.clear();
    smpl_full_loglike_.reserve(sample_size);
    smpl_full_loglike_.assign(sample_size, hsp::na::DOUBLE);

    // initialize  other diagnostic vectors
    // UPDATE: I need to store also power posterior sample
    // NOTE: it will depend on chosen rungs in that case
    if (!setup.is_power_post())
    {
        dataset.init_draws(sample_size);
        param.init_draws(sample_size);
        augm.init_draws(sample_size);
    }
}

/// Update Bayesian posterior sample.
void MCMC::Impl::update_posterior_sample_()
{
    // store sample of the full loglikelihood
    smpl_full_loglike_[smpl_count_] = full_loglike_;

    if (!setup.is_power_post())
    {
        param.update_draws(smpl_count_);
        dataset.update_draws(smpl_count_);
        augm.update_draws(smpl_count_);
    }

    ++smpl_count_;
}

/// Deep clear pre-existing posterior draws.
void MCMC::Impl::clear_posterior_sample_()
{
    smpl_count_ = 0;
    smpl_full_loglike_.clear();
    smpl_full_loglike_.shrink_to_fit();

    dataset.clear_draws();
    param.clear_draws();
    augm.clear_draws();
}
