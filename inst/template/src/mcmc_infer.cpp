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

#include "mcmc_impl.h"

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
void MCMC::Impl::_init_move_counts()
{
    for (const par_key &p : param_standard_infer_list + param_cache_infer_list)
        param_map.at(p).reset_accept_rate();
}

/* propose a valid move (retrun true), otherwise update refusal (return false) */
bool MCMC::Impl::propose_valid_param_move(const par_key &p)
{
    Variable &this_par = param_map.at(p);

    // set new value of the parameter and store the old to buffer
    this_par.propose(setup.prng());

    // recover buffer if not valid move
    if (!this_par.valid_move())
    {
        this_par.recover_value();
        return false;
    }

    return true;
}

/* Restore the parameter value in the buffer. */
void MCMC::Impl::recover_parameter(const par_key &p, bool is_cache_depen)
{
    // recover buffer value of parameter
    param_map.at(p).recover_value();

    // recover buffer value of full loglike
    full_loglike = full_buffer;

    // recover buffered loglikelihood in households
    for (Household &house : house_vector)
    {
        house.recover_buffer_values(is_cache_depen);
        if (is_augm_sample())
            for (Individual &member : house.infected())
                member.recover_augm_model_loglike();
    }

    recover_infec_profile(is_cache_depen);
}

/* Metropolis-Hastings to accept or not current parameter move. */
bool MCMC::Impl::accept_valid_param_move(const par_key &p)
{
    const double D_ll = full_loglike - full_buffer;
    const double move_logratio = param_map.at(p).get_move_logratio();
    return setup.unif() < std::exp(setup.inv_T * D_ll + move_logratio);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  MCMC::Impl : methods (diagnostics) ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init posterior sampling. UPDATE: */
void MCMC::Impl::_init_infer_diagnostics()
{
    ulong sample_size = setup.c_sample_size * setup.n_chains;

    // diagnostics over full loglikelihood
    smpl_count = 0;
    smpl_full_loglike.clear();
    smpl_full_loglike.reserve(sample_size);
    smpl_full_loglike.assign(sample_size, HSP_NA_DOUBLE);

    // initialize  other diagnostic vectors
    // UPDATE: I might want to store also power posterior sample
    if (!setup.is_power_post())
    {
        // sample over households
        smpl_house_loglike.clear();
        smpl_house_loglike.resize(restrict_hh_indexes.size());
        for (size_t R = 0; R < restrict_hh_indexes.size(); ++R)
        {
            smpl_house_loglike[R].clear();
            smpl_house_loglike[R].reserve(sample_size);
            smpl_house_loglike[R].assign(sample_size, HSP_NA_DOUBLE);
        }

        // sample over parameters
        for (const par_key &p : param_standard_infer_list + param_cache_infer_list)
        {
            smpl_param_value[p].clear();
            smpl_param_value[p].reserve(sample_size);
            smpl_param_value[p].assign(sample_size, HSP_NA_DOUBLE);
        }

        // sample over augmented variables
        if (!setup.id_trace.empty())
        {
            smpl_augm_value.clear();
            smpl_augm_value.resize(smpl_individual_star.size());
            for (size_t k = 0; k < smpl_individual_star.size(); ++k)
            {
                smpl_augm_value[k].clear();
                smpl_augm_value[k].reserve((augm_list_infer + augm_list_func).size());
                for (const augm_key &a : augm_list_infer + augm_list_func)
                {
                    smpl_augm_value[k][a].clear();
                    smpl_augm_value[k][a].reserve(sample_size);
                    smpl_augm_value[k][a].assign(sample_size, HSP_NA_DOUBLE);
                }
            }
        }
    }
}

/* Update Bayesian posterior sample. */
void MCMC::Impl::_update_posterior_sample()
{
    // store sample of the full loglikelihood
    smpl_full_loglike[smpl_count] = full_loglike;

    if (!setup.is_power_post())
    {
        // store sample for all parameters to be inferred
        for (const par_key &p : param_standard_infer_list + param_cache_infer_list)
            smpl_param_value.at(p)[smpl_count] = get_param(p);

        // store sample for restricted households
        for (size_t R = 0; R < restrict_hh_indexes.size(); ++R)
        {
            Household &house = house_vector[restrict_hh_indexes[R]];
            double x = house.get_obs_loglike() / house.multiplicity();
            for (Individual &member : house.infected())
                x += member.get_augm_model_loglike();
            smpl_house_loglike[R][smpl_count] = x;
        }

        // store sample for individual augmented variable followup
        if (!setup.id_trace.empty())
            for (size_t k = 0; k < smpl_individual_star.size(); ++k)
            {
                auto &smpl_augm_value_k = smpl_augm_value[k];
                Individual *member_star = smpl_individual_star[k];
                for (const augm_key &a : augm_list_infer + augm_list_func)
                    smpl_augm_value_k[a][smpl_count] = member_star->augm(a);
            }
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
void MCMC::Impl::_update_full_loglike(bool cache_depen)
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
            new_ll += _update_house_loglike(R, is_augm_ll, cache_depen);
        }
        catch (const std::exception &e)
        {
#pragma omp critical
            {
                if (!has_error.exchange(true))
                    error_message = "While computing loglikelihood for house " +
                                    house_vector[restrict_hh_indexes[R]].house_id() + ":\n" +
                                    std::string(e.what());
            }
        }
    }
    if (!error_message.empty())
        throw std::runtime_error(error_message);

    set_full_loglike(new_ll);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : update house loglikelihood  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Update and return contribution of given household to full logliklehood.
 * NOTE: house multiplicity doesn't work with data augmentation
 */
double MCMC::Impl::_update_house_loglike(size_t R, bool is_augm_ll, bool cache_depen)
{
    Household &house = house_vector[restrict_hh_indexes[R]];

    // buffer all household cache dependencies
    if (cache_depen)
    {
        house.buffer_all_cache();
        for (Individual &member : house.infected())
            _update_cache_dependecies(house, member);
    }

    // update household observation loglikelihood
    double output = _compute_and_set_household_obs_loglike(house);

    // update individual augmented loglikelihood
    if (is_augm_ll)
        for (Individual &member : house.infected())
            output += member.compute_and_set_augm_model_loglike();

    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : parameter move  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Attempt to move a single parameter for the given inverse temperature. */
bool MCMC::Impl::_move_single_parameter(const par_key &p, bool cache_depen)
{
    bool is_accepted = false;
    // proceed only if the new value is consistent
    // otherwise automatically execute refuse_move(par_key)
    if (propose_valid_param_move(p))
    {
        // if cache depending on parameter, update and buffer cache
        _update_full_loglike(cache_depen);

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
void MCMC::Impl::_move_all()
{
    // move augmented variables if any
    _move_all_augm();

    // move all simple model parameters
    for (const par_key &p : param_standard_infer_list)
        _move_single_parameter(p, false);

    // move parameters with cache dependecies
    for (const par_key &p : param_cache_infer_list)
        _move_single_parameter(p, true);
}

/* Move all augmented variables and model parameters and adapt their rule arguments. */
void MCMC::Impl::_move_all_and_adapt_params(ulong step)
{
    // move augmented variables if
    _move_all_augm();

    // move all model parameters and adapt their rule arguments
    const double ratio = static_cast<double>(step) / setup.n_adapt;

    // move all simple model parameters
    for (const par_key &p : param_standard_infer_list)
    {
        Variable &this_param = param_map.at(p);
        this_param.update_accept_rate(_move_single_parameter(p, false));
        this_param.adapt_rule(ratio);
    }

    // move parameters with cache dependecies
    for (const par_key &p : param_cache_infer_list)
    {
        Variable &this_param = param_map.at(p);
        this_param.update_accept_rate(_move_single_parameter(p, true));
        this_param.adapt_rule(ratio);
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
void MCMC::Impl::_run_monte_carlo()
{
    //==========================//
    //===== Initialization =====//
    //==========================//

    ulong &chain = setup.chain; // current chain
    _init_infer_diagnostics();  // initialize diagnostics

    // run multiple MCMC chains&
    for (chain = 1; chain <= setup.n_chains; ++chain)
    {
        // update progressbar text
        setup.prog_set_text();      // update progress bar text
        init_parameters();          // init parameter inference
        init_augmented_variables(); // init augmented dependencies
        _update_full_loglike(true); // init cache + full loglikelihood
        _init_move_counts();        // init move counts for adaptive

        // (Finite) Adaptative Iterations
        for (ulong i = 0; i < setup.n_adapt; i += setup.thin, setup.prog_update())
            for (ulong t = 1; t <= setup.thin; ++t)
                _move_all_and_adapt_params(i + t);

        // Burn-in Iterations
        for (ulong i = 0; i < setup.n_burnin; i += setup.thin, setup.prog_update())
            for (ulong t = 1; t <= setup.thin; ++t)
                _move_all();

        // Posterior sampling
        for (ulong i = setup.n_burnin; i < setup.n_steps; i += setup.thin, setup.prog_update())
        {
            for (ulong t = 1; t <= setup.thin; ++t)
                _move_all();
            _update_posterior_sample(); // store sample
        }
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC::Impl : evaluate loglikelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize model and compute the loglikelihood. */
double MCMC::Impl::evaluate()
{
    try
    {
        assert_proper_dataset();         // check if dataset is empty
        setup.n_chains = 0;              // use 0 chains as a state variable
        setup.rand_init_params = false;  // disable random initialization // UPDATE
        init_parameters();               // check no undefined parameters
        init_augmented_variables();      // initialize augmented dependencies
        infec_profile.assign_template(); // assign infectivity profile template
        _init_model_for_inference();     // initialize model for inference
        _update_full_loglike(true);      // initialize cache + full loglikelihood

        return full_loglike;
    }
    catch (const std::exception &e)
    {
        tools::error(e.what());
    }

    return HSP_NA_DOUBLE;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC::Impl : inference  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Performs the inference using Metropolis-Hastings. */
void MCMC::Impl::infer()
{

    try
    {
        setup.prog_declare();            // declare current progress bar
        assert_proper_dataset();         // check if dataset is empty
        assert_proper_mcmc_setup();      // check setup for mcmc sampling
        infec_profile.assign_template(); // assign infectivity profile template
        _init_model_for_inference();     // user initialization model

        setup.prog_open(setup.c_n_sweeps * setup.n_chains); // open progress bar
        _run_monte_carlo();                                 // run Monte Carlo
        setup.prog_close();                                 // close progress bar
    }
    catch (std::exception &e)
    {
        setup.prog_crash();
        tools::error(e.what());
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : power posteriors  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Power posteriors numerical integration for evidence.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
const std::vector<dbl_vec> MCMC::Impl::_power_posteriors(const std::string &scheme, ulong n_rungs)
{
    try
    {
        //==========================//
        //===== INITIALIZATION =====//
        //==========================//

        const str_vec options = {"adaptive", "power-fractions"};
        if (!contains(options, scheme))
            throw std::invalid_argument("Unknown power posterior scheme \"" + scheme + "\"." +
                                        tools::find_closest(scheme, options, 3));

        assert_proper_dataset();    // check if dataset is empty
        assert_proper_mcmc_setup(); // check setup for mcmc sampling

        // fix internal number of rungs (for progress bar)
        if (n_rungs >= 2)
            setup.n_rungs = n_rungs;
        else
            throw std::invalid_argument("N. of rungs must be greater than 1.");
        ulong &rung = setup.rung; // current rung

        // declare current progress bar
        setup.prog_declare();

        // initialize model for inference
        _init_model_for_inference();

        // implementation setup aliases
        double &t = setup.inv_T; // current inverse temperature

        //----- Buffer user setup -----//
        double buffer_inv_T = t;
        bool buffer_rand_init_params = setup.rand_init_params;

        // thermodynamic ladder estimates
        std::pair<double, double> stats = {}; // results of stats computation
        dbl_vec itm = {};                     // inverse temperature
        dbl_vec avg = {};                     // mean loglikelihood
        dbl_vec var = {};                     // variance loglikelihood

        //=================================//
        //===== INVERSE TEMPERATURE 0 =====//
        //=================================//

        setup.prog_open((setup.c_n_sweeps * setup.n_chains) * (setup.n_rungs + 1));

        // UPDATE: should I use previous inference for itm = 1 if available?
        // UPDATE: always under-estimate time if starting from itm = 0

        setup.rand_init_params = true;

        rung = 0;
        t = 0.0;
        _run_monte_carlo();
        stats = _smpl_full_loglike_stats();

        itm.push_back(t);
        avg.push_back(stats.first);
        var.push_back(stats.second);

        setup.rand_init_params = false;

        if (scheme == "adaptive")
        {
            //===========================//
            //===== ADAPTIVE LADDER =====//
            //===========================//

            //----- Inverse Temperature 1 -----//
            rung = 1;
            t = 1.0;
            _run_monte_carlo();
            stats = _smpl_full_loglike_stats();

            itm.push_back(t);
            avg.push_back(stats.first);
            var.push_back(stats.second);

            //----- Inverse Temperatures (0, 1) -----//
            size_t k = 0; // index to interval lower bound
            bool is_right = false;
            while (rung < n_rungs)
            {
                ++rung;

                //----- adaptive choice of inverse temperature -----//
                // intersection between tangents
                double new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
                // otherwise weighted average
                if ((new_t <= itm[k]) || (new_t >= itm[k + 1]))
                    new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);

                //----- power posterior sampling -----//
                t = new_t;
                _run_monte_carlo();
                stats = _smpl_full_loglike_stats();

                //----- find larger rectangular discretization error -----//
                is_right = (itm[k + 1] - t) * (avg[k + 1] - stats.first) > (new_t - itm[k]) * (stats.first - avg[k]);

                //----- insert in-between for all lists -----//
                itm.insert(itm.begin() + k + 1, t);
                avg.insert(avg.begin() + k + 1, stats.first);
                var.insert(var.begin() + k + 1, stats.second);

                // if right interval is chosen, index points to inserted as lower bound
                k += (is_right);
            }
        }
        else if (scheme == "power-fractions")
        {
            //==================================//
            //===== POWER FRACTIONS LADDER =====//
            //==================================//

            //----- Inverse Temperatures (0, 1] -----//
            for (rung = 1; rung <= n_rungs; ++rung)
            {
                // inverse tempearture = (rung / #rungs)^5
                t = std::pow(static_cast<double>(rung) / n_rungs, 5);
                _run_monte_carlo();
                stats = _smpl_full_loglike_stats();

                itm.push_back(t);
                avg.push_back(stats.first);
                var.push_back(stats.second);
            }
        }

        //----- Close progress bar -----//
        setup.prog_close();

        //----- Restore user setup -----//
        t = buffer_inv_T;
        setup.rand_init_params = buffer_rand_init_params;

        return {itm, avg, var};
    }
    catch (const std::exception &e)
    {
        setup.prog_crash();
        throw e;
    }
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : impl : simulation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Auxiliary method to initialize simulations. */
void MCMC::Impl::_aux_init_simulate()
{
    assert_proper_dataset(); // check if dataset is empty
    setup.n_chains = 0;      // use 0 chains as a state variable

    setup.n_threads = 1; // impose no parallelization
    setup.set_n_threads();

    setup.rand_init_params = false;  // disable random parameter initialization
    init_parameters();               // check no undefined parameters
    infec_profile.assign_template(); // assign infectivity profile template
    _init_model_for_simulation();    // initialize model for simulation
}

#ifdef _HOUSEPI_RCPP
/* Perform simulation and return simulated dataset. */
Rcpp::DataFrame MCMC::Impl::simulate()
{
    using namespace Rcpp;

    try
    {
        setup.prog_declare("simulate");

        const size_t n_cols = n_individuals;
        const str_vec _aux_augm_headers = subtr(augm_list_full,
                                                {AUX_IP_INST_PROFILE, AUX_IP_CUMUL_PROFILE});
        const str_vec _aux_dataset_headers = subtr(house_vector[0].get_member(0).keys_dataset(),
                                                   _aux_augm_headers);

        Rcpp::CharacterVector df_index = CharacterVector(n_cols);
        Rcpp::CharacterVector df_house_index = CharacterVector(n_cols);
        Rcpp::IntegerVector df_status = IntegerVector(n_cols);
        std::vector<Rcpp::NumericVector> df_dbl = {};

        // other data and augmented variable columns
        for (size_t s = 0; s < _aux_dataset_headers.size() + _aux_augm_headers.size(); ++s)
            df_dbl.emplace_back(n_cols);

        //----- initialization -----//
        _aux_init_simulate();
        setup.prog_open(n_houses());

        //----- run simulation -----//
        size_t row = 0;
        for (Household house : house_vector) // deep copy of the household
            try
            {
                // inplace simlation
                _simulate_single_household(house);

                // store simulation
                for (Individual &member : house)
                {
                    df_index[row] = member.id();
                    df_house_index[row] = house.house_id();
                    df_status[row] = member.status();

                    size_t col = 0;
                    for (const data_key &d : _aux_dataset_headers)
                    {
                        double x = member.data(d);
                        df_dbl[col++][row] = tools::is_na(x) ? NA_REAL : x;
                    }

                    for (const augm_key &a : _aux_augm_headers)
                    {
                        double x = member.augm(a);
                        df_dbl[col++][row] = tools::is_na(x) ? NA_REAL : x;
                    }

                    ++row;
                }
                setup.prog_update(); // update progress bar
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error("While simulating house " +
                                         house.house_id() + " :\n" +
                                         std::string(e.what()));
            }

        // close progress bar
        setup.prog_close();

        // assign vectors to dataframe
        DataFrame output = DataFrame::create(
            _["id"] = std::move(df_index),
            _["house_id"] = std::move(df_house_index),
            _["status"] = std::move(df_status));

        size_t col = 0;
        for (const data_key &d : _aux_dataset_headers)
            output.push_back(std::move(df_dbl[col++]), d);

        for (const augm_key &a : _aux_augm_headers)
            output.push_back(std::move(df_dbl[col++]), a);

        return output;
    }
    catch (std::exception &e)
    {
        setup.prog_crash();
        tools::error(e.what());
    }

    return DataFrame::create();
}

Rcpp::List MCMC::Impl::estimate_evidence(const std::string &scheme, ulong n_rungs)
{
    using namespace Rcpp;

    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = _power_posteriors(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (ulong k = 0; k < n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return List::create(_["estimate"] = output,
                            _["scheme"] = scheme,
                            _["n_rungs"] = n_rungs,
                            _["inv_temper"] = wrap(inv_tem),
                            _["mean_loglike"] = wrap(avg_vec),
                            _["var_loglike"] = wrap(var_vec));
    }
    catch (const std::exception &e)
    {
        tools::error("While estimating evidence with power posteriors:\n" +
                     std::string(e.what()));
    }
    return List();
}

#else

/* Print the simulated dataset to the provided file. */
void MCMC::Impl::print_simulated_dataset(const std::string &pathtofile)
{

    try
    {
        //----- initialization -----//
        setup.prog_declare("simulate");
        _aux_print_dataset_init_(pathtofile);
        _aux_init_simulate();

        setup.prog_open(n_houses()); // open progress bar

        //----- run simulation -----//
        // note: deep copy of the household
        for (Household house : house_vector)
            try
            {
                _simulate_single_household(house); // inplace simulation
                _aux_print_dataset_update_(house); // print the simulated dataset
                setup.prog_update();               // update progress bar
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While simulating house \"" + house.house_id() + "\" :\n" + e.what());
            }

        setup.prog_close();             // close progress bar
        _aux_print_dataset_finalize_(); // finalize output
    }
    catch (std::exception &e)
    {
        setup.prog_crash();
        tools::error(e.what());
    }
}

/* Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
    -----
    Ref: Friel, N., Hurn, M. & Wyse, J.
    Improving power posterior estimation of statistical evidence.
    Stat Comput 24, 709–723 (2014).
*/
void MCMC::Impl::print_evidence(const std::string &scheme, ulong n_rungs)
{
    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = _power_posteriors(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (ulong k = 0; k < n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        // print result to stdout
        tools::print_stdout(std::to_string(output) + "\n");
    }
    catch (const std::exception &e)
    {
        tools::error("While estimating evidence with power posteriors:\n" +
                     std::string(e.what()));
    }
}

#endif
