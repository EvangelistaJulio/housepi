/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_build_param.cpp
 *
 */

#include <atomic>

#include "mcmc_build.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Assign parameter definitions to dictionary. */
void MCMC::Impl::build_param()
{
    try
    {
        // reference to parameter definition table
        const dict_str_vec &deftable = _param_deftable;

        if (!deftable.empty())
        {
            // check entries
            str_vec headers = tools::headers_from_deftable(deftable);
            tools::check_for_missing_entries(headers, param_file_headers);
            tools::check_for_unused_entries(headers, param_file_headers);

            // drop current list of parameters
            param_map.clear();
            param_const_list.clear();
            param_standard_infer_list.clear();
            param_cache_infer_list.clear();

            //===== Load parameters from the row with options =====//
            for (size_t k = 0; k < deftable.at("name").size(); ++k)
            {
                const std::string &name = deftable.at("name")[k]; // parameter name
                double value = HOUSEPI_NA_DOUBLE;                 // init value
                parID p = name;

                // Add entry to shared parameter dictionary -----//
                if (param_map.count(p))
                {
                    // parameter is already defined
                    throw std::invalid_argument(
                        "Parameter \"" + p + "\" has multiple definitions.");
                }
                else if (tools::stobool(deftable.at("const")[k]))
                {
                    // assign constant parameter to the set of variables
                    if (tools::is_na(deftable.at("value")[k]))
                        throw std::invalid_argument(
                            "Undefined value for constant parameter \"" + name + "\".");
                    else
                        value = tools::stod(deftable.at("value")[k], "parameter");

                    // assign constant parameter to the set of variables
                    param_map.emplace(p, Variable(name, value));
                    // update the lists of constant parameters
                    param_const_list.push_back(p);
                }
                else
                {
                    // load value
                    if (!tools::is_na(deftable.at("value")[k]))
                        value = tools::stod(deftable.at("value")[k], "parameter");

                    // load boolean flag for discrete variables
                    bool is_discrete = false;
                    if (!tools::is_na(deftable.at("discr")[k]))
                        is_discrete = tools::stobool(deftable.at("discr")[k]);

                    // load prior name
                    const std::string &prior_name = deftable.at("prior")[k];

                    // load rule name
                    std::string rule_name = "draw";
                    if (!tools::is_na(deftable.at("rule")[k]))
                        rule_name = deftable.at("rule")[k];

                    // load rule scale
                    double rule_scale = 1.;
                    if (!tools::is_na(deftable.at("scale")[k]))
                        rule_scale = tools::stod(deftable.at("scale")[k]);

                    //----- Create a new Parameter Variable -----//
                    param_map.emplace(p, Variable(name, value, is_discrete));
                    Variable &new_param = param_map.at(p);

                    // assign a prior distribution to the parameter
                    new_param.assign_distr(prior_name);

                    // assign scalar arguments to prior distribution
                    str_vec arg_key = tools::parse_list_args(deftable.at("arg_key")[k]);
                    str_vec arg_val = tools::parse_list_args(deftable.at("arg_val")[k]);
                    for (ushort c = 0; c < arg_key.size(); ++c)
                        // UPDATE: is this check really necessary?
                        if (!tools::is_na(arg_key[c]))
                        {
                            const std::string &x = arg_val[c];
                            new_param.assign_distr_attr(
                                arg_key[c], ATTR_ARGUMENT,
                                (tools::is_na(x) ? HOUSEPI_NA_DOUBLE : tools::stod(x)));
                        }

                    // load parameter range
                    double lower = HOUSEPI_NA_DOUBLE, upper = HOUSEPI_NA_DOUBLE;
                    if (!tools::is_na(deftable.at("min")[k]))
                        lower = tools::stod(deftable.at("min")[k], "parameter minimum");
                    if (!tools::is_na(deftable.at("max")[k]))
                        upper = tools::stod(deftable.at("max")[k], "parameter maximum");

                    // set the range of values to the parameter
                    if (!tools::is_na(lower))
                        new_param.assign_distr_attr("min", ATTR_RANGE, lower);
                    if (!tools::is_na(upper))
                        new_param.assign_distr_attr("max", ATTR_RANGE, upper);
                    new_param.init_distr_reshaping();

                    // assign the rule for the proposer to the parameter
                    new_param.assign_proposer(rule_name, rule_scale);

                    // update the lists of parameters to be inferred
                    param_standard_infer_list.push_back(p);
                }
            }
        }
    }
    catch (const std::exception &e)
    {
        tools::error("While building parameters:\n" +
                     std::string(e.what()));
    }
}

/* Randomly initialize values of the parameters to be inferred */
void MCMC::Impl::init_parameters()
{
    for (const parID &p : param_standard_infer_list + param_cache_infer_list)
    {
        Variable &this_par = param_map.at(p);
        // user initialization (for multiple chains)
        this_par.user_reset();

        if (setup_rand_init_params)
        {
            if (tools::is_na(this_par.show_distr_name()))
                throw std::invalid_argument(
                    "Parameters with undefined values are not allowed.");
            else
                this_par.rand(_prng);
        }
        else if (tools::is_na(this_par.get()))
            throw std::invalid_argument(
                "Value of parameter \"" + p + "\" is undefined with no random initialization.");
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC::Impl : methods (parameter moves)  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set accept rates to zero. */
void MCMC::Impl::init_move_counts()
{
    for (const parID &p : param_standard_infer_list + param_cache_infer_list)
        param_map.at(p).reset_accept_rate();
}

/* propose a valid move (retrun true), otherwise update refusal (retrun false) */
bool MCMC::Impl::propose_valid_param_move(const parID &p)
{
    Variable &this_par = param_map.at(p);

    // set new value of the parameter and store the old to buffer
    this_par.propose(_prng);

    // recover buffer if not valid move
    if (!this_par.valid_move())
    {
        this_par.recover_value();
        return false;
    }

    return true;
}

/* Restore the parameter value in the buffer. */
void MCMC::Impl::recover_parameter(const parID &p, bool cache_depen)
{
    // recover buffer value of parameter
    param_map.at(p).recover_value();

    // recover buffer value of full loglike
    full_loglike = full_buffer;

    // recover buffered loglikelihood in households
    for (Household &house : house_vector)
    {
        house.recover_loglike_value();
        if (cache_depen)
            house.recover_all_cache();
    }

    if (cache_depen)
        recover_infec_profile();
}

/* Update the acceptance rate of the paramter. */
void MCMC::Impl::update_param_accept_rate(const parID &p, bool is_accepted)
{
    param_map.at(p).update_accept_rate(is_accepted);
}

/* Adapt the parameter Metropolis move. */
void MCMC::Impl::adapt_param_rule(const parID &p, double ratio)
{
    param_map.at(p).adapt_rule(ratio);
}

/* Metropolis-Hastings to accept or not current parameter move. */
bool MCMC::Impl::accept_valid_param_move(const parID &p)
{
    const double D_ll = full_loglike - full_buffer;
    const double move_logratio = param_map.at(p).get_move_logratio();
    return uniform() < std::exp(setup_inv_T * D_ll + move_logratio);
}

/* Getter for full loglikelihood. */
double MCMC::Impl::get_full_loglike() const
{
    return full_loglike;
}

/* Set new full loglikelihood value and buffer old value. */
void MCMC::Impl::set_full_loglike(double new_ll)
{
    full_buffer = full_loglike;
    full_loglike = new_ll;
}

/* Increment full loglikelohood by a delta, buffer old value. */
void MCMC::Impl::increment_full_loglike(double delta)
{
    set_full_loglike(full_loglike + delta);
}

/* Compute statistics from sample of full loglikelihood. */
const std::pair<double, double> MCMC::Impl::smpl_full_loglike_stats() const
{
    return tools::compute_stats(smpl_full_loglike);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  MCMC::Impl : methods (diagnostics) ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init posterior sampling. */
void MCMC::Impl::init_diagnostics()
{
    ulong sample_size = setup_c_sample_size * setup_n_chains;

    // diagnostics over full loglikelihood
    smpl_count = 0;
    smpl_full_loglike.clear();
    smpl_full_loglike = dbl_vec(sample_size);

    // initialize  other diagnostic vectors
    // UPDATE: I might want to store also power posterior sample
    if (!is_power_post())
    {
        // diagnostics over parameters
        for (const parID &p : param_standard_infer_list + param_cache_infer_list)
            smpl_param_value[p] = dbl_vec(sample_size, HOUSEPI_NA_DOUBLE);

        // diagnostics over households
        // UPDATE: init a=only for selected households
        for (size_t H : restrict_hh_indexes)
            house_vector[H].init_diagnostics(sample_size);

        // diagnostics over augmented variables
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_list_infer + augm_list_func)
                smpl_augm_value[a] = dbl_vec(sample_size, HOUSEPI_NA_DOUBLE);
    }
}

/* Update Bayesian diagnostics */
void MCMC::Impl::update_diagnostics()
{
    // diagnostics for full loglike
    smpl_full_loglike[smpl_count] = get_full_loglike();

    if (!is_power_post())
    {
        // diagnostics for all parameters
        for (const parID &p : param_standard_infer_list + param_cache_infer_list)
            smpl_param_value.at(p)[smpl_count] = get_param(p);

        // diagnostics for selected households
        for (size_t H : restrict_hh_indexes)
            house_vector[H].update_diagnostics(smpl_count);

        // diagnostics for individual followup
        if (!tools::is_na(setup_id_trace))
            for (const augmID &a : augm_list_infer + augm_list_func)
                smpl_augm_value.at(a)[smpl_count] = smpl_individual_star->augm(a);
    }

    ++smpl_count;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : update house loglikelihood  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute new loglikelihood, store old value in buffer. Same for households. */
void MCMC::Impl::update_full_loglike(const bool cache_depen)
{
    const bool is_augm_ll = is_augm_model();
    double new_ll = 0.0;
    std::atomic<bool> has_error(false);

    // buffer infectivity profile if lookup tables are active
    if (cache_depen)
        update_infec_profile();

#pragma omp parallel for reduction(+ : new_ll)
    for (size_t R = 0; R < restrict_hh_indexes.size(); ++R)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        try
        {
            new_ll += update_house_loglike(R, is_augm_ll, cache_depen);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    _error_message = "While computing loglikelihood for house " +
                                     house_vector[restrict_hh_indexes[R]].house_id() + ":\n" +
                                     std::string(e.what());
            }
        }
    }
    if (!_error_message.empty())
        tools::error(_error_message);

    set_full_loglike(new_ll);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : update house loglikelihood  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Update and return contribution of given household to full logliklehood. */
double MCMC::Impl::update_house_loglike(size_t R, const bool is_augm_ll, const bool cache_depen)
{
    Household &house = house_vector[restrict_hh_indexes[R]];

    // buffer all household cache dependencies
    if (cache_depen)
    {
        house.buffer_all_cache();
        for (const indiv &i : house.get_infected_status_hhindexes())
            _update_cache_dependecies(house, i);
    }

    // update household observation loglikelihood
    double output = _compute_house_obs_loglike(house);
    house.set_obs_loglike(output);

    // update household augmented loglikelihood
    // NOTE: house multiplicity doesn't work with data augmentation
    if (is_augm_ll)
    {
        house.set_augm_loglike(_compute_augm_model_loglike(house));
        output += house.get_augm_loglike();
    }

    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : parameter move  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Attempt to move a single parameter for the given inverse temperature. */
bool MCMC::Impl::move_single_parameter(const parID &p, bool cache_depen)
{
    bool is_accepted = false;
    // proceed only if the new value is consistent
    // otherwise automatically execute refuse_move(parID)
    if (propose_valid_param_move(p))
    {
        // if cache depending on parameter, update and buffer cache
        update_full_loglike(cache_depen);

        // perform Metropolis-Hastings to accept or not move
        if (accept_valid_param_move(p))
            is_accepted = true;
        else
            recover_parameter(p, cache_depen);
    }

    return is_accepted;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Move all model parameters. */
void MCMC::Impl::move_all()
{
    // move augmented variables if any
    move_all_augm();

    // move all simple model parameters
    for (const parID &p : param_standard_infer_list)
        move_single_parameter(p, false);

    // move parameters with cache dependecies
    for (const parID &p : param_cache_infer_list)
        move_single_parameter(p, true);
}

/* Move all augmented variables and model parameters and adapt their rule arguments. */
void MCMC::Impl::move_all_and_adapt_parameters(ulong step)
{
    // move augmented variables if
    move_all_augm();

    // move all model parameters and adapt their rule arguments
    double adapt_ratio = static_cast<double>(step) / setup_n_adapt;

    // move all simple model parameters
    for (const parID &p : param_standard_infer_list)
    {
        update_param_accept_rate(p, move_single_parameter(p, false));
        adapt_param_rule(p, adapt_ratio);
    }

    // move parameters with cache dependecies
    for (const parID &p : param_cache_infer_list)
    {
        update_param_accept_rate(p, move_single_parameter(p, true));
        adapt_param_rule(p, adapt_ratio);
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : run markov chain  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Markov chain core, used for inference and power posteriors.
    note: sweep is the number of previously completed thin sweeps. */
void MCMC::Impl::run_monte_carlo(progress::Bar &progressbar)
{
    //==========================//
    //===== Initialization =====//
    //==========================//

    init_diagnostics(); // initialize diagnostics

    // run multiple MCMC chains&
    for (chain = 1; chain <= setup_n_chains; ++chain)
    {
        // update progressbar text
        progressbar.set_text(fmt_progtxt());
        init_parameters();          // init parameter inference
        init_augmented_variables(); // init augmented dependencies
        update_full_loglike(true);  // init cache + full loglikelihood
        init_move_counts();         // init move counts for adaptive

        // (Finite) Adaptative Iterations
        for (ulong i = 0; i < setup_n_adapt; i += setup_thin, progressbar.update())
            for (ulong t = 1; t <= setup_thin; ++t)
                move_all_and_adapt_parameters(i + t);

        // Burn-in Iterations
        for (ulong i = 0; i < setup_n_burnin; i += setup_thin, progressbar.update())
            for (ulong t = 1; t <= setup_thin; ++t)
                move_all();

        // Posterior sampling
        for (ulong i = setup_n_burnin; i < setup_n_steps; i += setup_thin, progressbar.update())
        {
            for (ulong t = 1; t <= setup_thin; ++t)
                move_all();
            update_diagnostics(); // store sample
        }
    }
}
