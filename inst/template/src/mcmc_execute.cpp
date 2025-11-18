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
 * mcmc_execute.cpp
 *
 */

#include "tools/string.h"
#include "tools/error.h"
#include "mcmc_impl.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC::Impl : evaluate log-likelihood  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize model and evaluate log-likelihood.
double MCMC::evaluate_log_lik(const std::string &house_id)
{
    try
    {
        pImpl_->init_evaluate();
        // mode for computing log-likelihood of specific household
        if (!hsp::na::is(house_id))
            return pImpl_->evaluate_single_house_log_lik(house_id);
        // initialize cached values and then compute total log-likelihood
        return pImpl_->update_total_log_lik();
    }
    catch (const std::exception &e)
    {
        hsp::msg::error(e.what());
    }
}

/// Assert setup definitions are ok to proceed with log-likelihood evaluation.
void MCMC::Impl::init_evaluate()
{
    setup.n_chains = 0;            // use 0 chains as a state variable
    assert_non_empty_data_wide_(); // check if data is empty
    warn_if_all_missing_status_(); // warn if all missing status
    param.user_reset();            // initialize non constant parameters
    augm.init();                   // initialize individual augmented relations
    augm.user_reset();             // initialize non constant augmented
    init_model_for_inference();    // initialize model for inference
}

/// Evaluate log-likelihood for the given household.
double MCMC::Impl::evaluate_single_house_log_lik(const std::string &house_id)
{
    atomic_error.reset();   // reset atomic error
    infec_profile.update(); // eventually buffer infectivity profile
    double output = atomic_materialize_house_log_lik_(data_wide(house_id), augm.is_any_model());
    atomic_error.throw_exception(); // catch and throw if previous failed
    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  MCMC : update house log-likelihood  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute new log-likelihood, store old value in buffer. Same for households.
double MCMC::Impl::update_total_log_lik(bool is_cache_depen)
{
    const bool is_augm_ll = augm.is_any_model();
    double new_ll = 0.0;

    // buffer infectivity profile if lookup tables are active
    if (is_cache_depen)
        infec_profile.update();

    atomic_error.reset();
#pragma omp parallel for reduction(+ : new_ll)
    for (size_t R = 0; R < multiplicity.indexes.size(); ++R)
    {
        if (atomic_error.is_pending())
            continue;

        new_ll += atomic_materialize_house_log_lik_(data_wide.restricted(R), is_augm_ll, is_cache_depen);
    }
    atomic_error.throw_exception();

    return new_ll;
}

/// Update and return contribution of given household to full logliklehood.
/// NOTE: house multiplicity doesn't work with data augmentation
double MCMC::Impl::atomic_materialize_house_log_lik_(Household &house, bool is_augm_ll, bool is_cache_depen)
{
    double output = 0;
    try
    {
        // buffer all household cache dependencies
        if (is_cache_depen)
        {
            house.buffer_all_cache();
            for (Individual &member : house.infected())
                update_cache_dependecies_(house, member);
        }

        // compute and set household observation log-likelihood
        output = materialize_house_obs_log_lik_(house);
        // compute and set individual augmented log-likelihood
        if (is_augm_ll)
            for (Individual &member : house.infected())
                output += member.materialize_augm_model_log_lik();
    }
    catch (const std::exception &e)
    {
        atomic_error.set_if_first("While computing log-likelihood for " +
                                  house.err_name() + ":\n" + std::string(e.what()));
    }
    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : Impl : inference  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Performs the inference using Metropolis-Hastings.
void MCMC::infer()
{
    Impl::Setup &setup = pImpl_->setup;
    try
    {
        setup.prog_declare();                               // declare current progress bar
        pImpl_->init_monte_carlo();                         // initialize mcmc sampling
        pImpl_->init_model_for_inference();                 // user initialization model
        pImpl_->augm.init();                                // initialize augmented variables
        setup.prog_open(setup.c_n_sweeps * setup.n_chains); // open progress bar

        pImpl_->run_monte_carlo(); // run Monte Carlo
        setup.prog_close();        // close progress bar
    }
    catch (std::exception &e)
    {
        setup.prog_crash();
        hsp::msg::error(e.what());
    }
}

/// Assert setup definitions are ok to proceed with mcmc sampling.
void MCMC::Impl::init_monte_carlo()
{
    setup.assert_mcmc_setup();     // assert if setup is ok for mcmc
    assert_non_empty_data_wide_(); // check if data is empty
    warn_if_all_missing_status_(); // warn if all missing status

    // check if there is actually something to sample
    if (!param.is_any_sample() && !augm.is_any_sample())
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
void MCMC::Impl::run_monte_carlo()
{
    unsigned long &chain = setup.chain; // current chain
    init_posterior_sample_();           // initialize diagnostics

    // run multiple MCMC chains&
    for (chain = 1; chain <= setup.n_chains; ++chain)
    {
        setup.prog_set_text(); // set progress bar text
        param.user_reset();    // reset parameter user definitions
        augm.user_reset();     // reset augmented user definitions
        set_total_log_lik_(update_total_log_lik()); // init cached, set total log-likelihood
        init_move_counts_(); // init move counts for adaptive

        // Finite adaptative iterations
        for (unsigned long i = 0; i < setup.n_adapt; i += setup.n_thin, setup.prog_update())
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
            {
                move_all_();
                adapt_proposal_of_sampled_params_(i + t);
            }

        // Burn-in iterations
        for (unsigned long i = 0; i < setup.n_burnin; i += setup.n_thin, setup.prog_update())
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_();

        // Posterior sampling
        for (unsigned long i = setup.n_burnin; i < setup.n_steps; i += setup.n_thin, setup.prog_update())
        {
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_();
            store_posterior_sample_();
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
    // move augmented variables to sample, if any
    if (augm.is_any_sample())
        move_all_augm_();

    // move all standard model parameters
    for (const std::string &p : param.list_adaptable_sample_standard)
        param.at(p).update_accept_rate(move_single_parameter_(p, false));
    for (const std::string &p : param.list_non_adapt_sample_standard)
        move_single_parameter_(p, false);

    // move parameters with cached dependecies
    for (const std::string &p : param.list_adaptable_sample_cached)
        param.at(p).update_accept_rate(move_single_parameter_(p, true));
    for (const std::string &p : param.list_non_adapt_sample_cached)
        move_single_parameter_(p, true);
}

/// Attempt to move a single parameter for the given inverse temperature.
bool MCMC::Impl::move_single_parameter_(const std::string &p, bool is_cache_depen)
{
    // proceed only if the new value is consistent
    // otherwise automatically refuse move
    if (propose_valid_param_move_(p))
    {
        // if parameter is cache dependent, update cache
        set_total_log_lik_(update_total_log_lik(is_cache_depen));

        // perform Metropolis-Hastings to accept or not move
        if (accept_valid_param_move_(p))
            return true; // accepted move
        recover_parameter_(p, is_cache_depen);
    }
    return false; // refused move
}

/// Adpat the proposal of model parameters.
void MCMC::Impl::adapt_proposal_of_sampled_params_(unsigned long step)
{
    const double ratio = static_cast<double>(step) / setup.n_adapt;

    // adpat all standard model parameters
    for (const std::string &p : param.list_adaptable_sample_standard)
        param.at(p).adapt_proposal(ratio);

    // adapt all parameters with cached dependecies
    for (const std::string &p : param.list_adaptable_sample_cached)
        param.at(p).adapt_proposal(ratio);
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

    // recover buffer value of full log-like
    total_log_lik_ = total_buffer_;

    // recover buffered log-likelihood in households
    for (Household &house : data_wide.restricted())
    {
        house.recover_buffer_values(is_cache_depen);
        if (augm.is_any_sample())
            for (Individual &member : house.infected())
                member.recover_augm_model_log_lik();
    }

    if (is_cache_depen)
        infec_profile.recover();
}

/// Metropolis-Hastings to accept or not current parameter move.
bool MCMC::Impl::accept_valid_param_move_(const std::string &p)
{
    const double D_ll = total_log_lik_ - total_buffer_;
    const double move_logratio = param.at(p).get_move_logratio();
    return setup.unif() < std::exp(setup.inv_T * D_ll + move_logratio);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : Impl : posterior sample  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize posterior sample of parameters and respective priors.
void MCMC::Impl::Parameter::init_draws(unsigned long sample_size)
{
    smpl_values.clear();
    smpl_values.reserve(list_sample.size());
    smpl_priors.clear();
    smpl_priors.reserve(list_sample.size());
    for (const std::string &p : list_sample)
    {
        smpl_values[p].clear();
        smpl_values[p].reserve(sample_size);
        smpl_values[p].assign(sample_size, hsp::na::DOUBLE);
        smpl_priors[p].clear();
        smpl_priors[p].reserve(sample_size);
        smpl_priors[p].assign(sample_size, hsp::na::DOUBLE);
    }
}

/// Store current sample of parameters and respective priors.
void MCMC::Impl::Parameter::store_draws(unsigned long smpl_count)
{
    for (const std::string &p : list_sample)
    {
        const Variable &this_par = this->at(p);
        smpl_values[p][smpl_count] = this_par.get();
        smpl_priors[p][smpl_count] = this_par.get_distr_logpdf();
    }
}

/// Deep clear existing posterior sample of parameters and respective priors.
void MCMC::Impl::Parameter::clear_draws()
{
    for (const std::string &p : list_sample)
    {
        smpl_values[p].clear();
        smpl_values[p].shrink_to_fit();
        smpl_priors[p].clear();
        smpl_priors[p].shrink_to_fit();
    }
    smpl_values.clear();
    smpl_values.rehash(0);
    smpl_priors.clear();
    smpl_priors.rehash(0);
}

/// Initialize posterior sample of parameters.
void MCMC::Impl::DataWide::init_draws(unsigned long sample_size)
{
    // sample over households
    smpl_values.clear();
    smpl_values.resize(pImpl_.multiplicity.indexes.size());
    for (size_t R = 0; R < pImpl_.multiplicity.indexes.size(); ++R)
    {
        smpl_values[R].clear();
        smpl_values[R].reserve(sample_size);
        smpl_values[R].assign(sample_size, hsp::na::DOUBLE);
    }
}

/// Store current sample of log-likelihood for restricted households.
void MCMC::Impl::DataWide::store_draws(unsigned long smpl_count)
{
    // store sample for restricted households
    for (size_t R = 0; R < pImpl_.multiplicity.indexes.size(); ++R)
    {
        Household &house = house_vector_[pImpl_.multiplicity.indexes[R]];
        double x = house.get_obs_log_lik() / house.multiplicity();
        for (const Individual &member : house.infected())
            x += member.get_augm_model_log_lik();
        smpl_values[R][smpl_count] = x;
    }
}

/// Deep clear posterior sample of pointwise log-likelihood.
void MCMC::Impl::DataWide::clear_draws()
{
    for (auto &vec : smpl_values)
    {
        vec.clear();
        vec.shrink_to_fit();
    }
    smpl_values.clear();
    smpl_values.shrink_to_fit();
}

/// Init posterior sampling.
void MCMC::Impl::init_posterior_sample_()
{
    unsigned long sample_size = setup.c_sample_size * setup.n_chains;

    // diagnostics over full log-likelihood
    smpl_count = 0;
    smpl_total_log_lik.clear();
    smpl_total_log_lik.reserve(sample_size);
    smpl_total_log_lik.assign(sample_size, hsp::na::DOUBLE);

    // initialize  other diagnostic vectors
    // UPDATE: I need to store also power posterior sample
    // NOTE: it will depend on chosen rungs in that case
    if (!setup.is_power_post())
    {
        data_wide.init_draws(sample_size);
        param.init_draws(sample_size);
        augm.init_draws(sample_size);
    }
}

/// Store current Bayesian posterior sample.
void MCMC::Impl::store_posterior_sample_()
{
    // store sample of the full log-likelihood
    smpl_total_log_lik[smpl_count] = total_log_lik_;

    if (!setup.is_power_post())
    {
        data_wide.store_draws(smpl_count);
        param.store_draws(smpl_count);
        augm.store_draws(smpl_count);
    }

    ++smpl_count;
}

/// Deep clear pre-existing posterior draws.
void MCMC::Impl::clear_posterior_sample_()
{
    smpl_count = 0;
    smpl_total_log_lik.clear();
    smpl_total_log_lik.shrink_to_fit();

    data_wide.clear_draws();
    param.clear_draws();
    augm.clear_draws();
}
void MCMC::clear_memory() { pImpl_->clear_memory(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : estimate log_evidence  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Welford online's algorithm for mean and variance.
/// Ref :
/// ----
/// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
std::pair<double, double> compute_stats(const dbl_vec &data)
{
    size_t i = 0;
    size_t N = data.size();
    double x, delta, delta2;
    double sqr_diff = 0.;

    double avg = 0.;
    while (i < N)
    {
        x = data[i++];
        delta = x - avg;
        avg += delta / i;
        delta2 = x - avg;
        sqr_diff += delta * delta2;
    }

    return {avg, sqr_diff / (N - 1.)};
}

/// Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
///
/// Ref
/// -----
/// Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
std::vector<dbl_vec> MCMC::Impl::power_posteriors(const std::string &scheme, unsigned long n_rungs)
{
    // check scheme choice
    const str_vec options = {"adaptive", "power-fractions"};
    if (!(options % scheme))
        throw std::invalid_argument(
            "Invoked invalid power posterior scheme '" + scheme + "'." +
            hsp::str::find_closest(scheme, options));

    // check number of rungs for the thermodynamic ladder
    if (n_rungs < 2)
        throw std::invalid_argument(
            "N. of rungs was not greater than 1.");
    setup.n_rungs = n_rungs; // NOTE: set before progress bar

    // initialize mcmc sampling
    init_monte_carlo();

    // buffer user setup
    double buffer_inv_T = setup.inv_T;

    setup.prog_declare();       // declare current progress bar
    init_model_for_inference(); // initialize model for inference
    augm.init();                // initialize augmented variables
    setup.prog_open(setup.c_n_sweeps * setup.n_chains * (setup.n_rungs + 1));

    // thermodynamic ladder estimates
    dbl_vec itm = {}; // inverse temperature
    dbl_vec avg = {}; // mean log-likelihood
    dbl_vec var = {}; // variance log-likelihood

    // run power posteriors
    power_post_zero_(itm, avg, var);
    if (scheme == "adaptive")
        power_post_adaptive_(itm, avg, var);
    else if (scheme == "power-fractions")
        power_post_fractions_(itm, avg, var);

    // close progress bar
    setup.prog_close();

    // restore user setup
    setup.inv_T = buffer_inv_T;

    return {itm, avg, var};
}

/// Power posteriors at 0 inverse temperature.
void MCMC::Impl::power_post_zero_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var)
{
    // UPDATE: should I use previous inference for itm = 1 if available?
    // UPDATE: always under-estimate time if starting from itm = 0

    setup.rung = 0;
    setup.inv_T = 0.0;
    run_monte_carlo();
    std::pair<double, double> stats = compute_stats(smpl_total_log_lik);

    itm.push_back(setup.inv_T);
    avg.push_back(stats.first);
    var.push_back(stats.second);
}

/// Power posteriors with adaptive scheme for the thermodynamic ladder.
void MCMC::Impl::power_post_adaptive_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var)
{
    // Inverse Temperature 1
    setup.rung = 1;
    setup.inv_T = 1.0;
    run_monte_carlo();
    std::pair<double, double> stats = compute_stats(smpl_total_log_lik);

    itm.push_back(setup.inv_T);
    avg.push_back(stats.first);
    var.push_back(stats.second);

    // Inverse Temperatures (0, 1)
    size_t k = 0; // index to interval lower bound
    bool is_right = false;
    while (setup.rung < setup.n_rungs)
    {
        ++setup.rung;

        // adaptive choice of inverse temperature
        // intersection between tangents
        double new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
        // otherwise weighted average
        if (new_t <= itm[k] || new_t >= itm[k + 1])
            new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);

        // power posterior sampling
        setup.inv_T = new_t;
        run_monte_carlo();
        stats = compute_stats(smpl_total_log_lik);

        // find larger rectangular discretization error
        is_right = (itm[k + 1] - setup.inv_T) * (avg[k + 1] - stats.first) > (new_t - itm[k]) * (stats.first - avg[k]);

        // insert in-between for all lists
        itm.insert(itm.begin() + k + 1, setup.inv_T);
        avg.insert(avg.begin() + k + 1, stats.first);
        var.insert(var.begin() + k + 1, stats.second);

        // if right interval is chosen, index points to inserted as lower bound
        k += (is_right);
    }
}

/// Power posteriors with power fractions scheme for the thermodynamic ladder.
void MCMC::Impl::power_post_fractions_(dbl_vec &itm, dbl_vec &avg, dbl_vec &var)
{
    for (setup.rung = 1; setup.rung <= setup.n_rungs; ++setup.rung)
    {
        // inverse tempearture = (rung / #rungs)^5
        setup.inv_T = std::pow(static_cast<double>(setup.rung) / setup.n_rungs, 5);
        run_monte_carlo();

        std::pair<double, double> stats = compute_stats(smpl_total_log_lik);
        itm.push_back(setup.inv_T);
        avg.push_back(stats.first);
        var.push_back(stats.second);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : Impl : standalone  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _HOUSEPI_RCPP

#include <fstream>
#include <iomanip>
#include "tools/read.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : Impl : read  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#define MAKE_READ_METHOD(module, section_name, csv_object)                \
    void MCMC::read_##module(const std::string &f)                        \
    {                                                                     \
        try                                                               \
        {                                                                 \
            pImpl_->module.build(hsp::read::csv_object(f));               \
        }                                                                 \
        catch (const std::exception &e)                                   \
        {                                                                 \
            hsp::msg::error(std::string("While reading ") +               \
                            section_name + " definitions:\n" + e.what()); \
        }                                                                 \
    }

MAKE_READ_METHOD(setup, "setup", csv_list)
MAKE_READ_METHOD(infec_profile, "infectivity profile", csv_list)
MAKE_READ_METHOD(param, "parameters", csv_table)
MAKE_READ_METHOD(data_wide, "data-wide", csv_table)
MAKE_READ_METHOD(multiplicity, "multiplicity", csv_table)
MAKE_READ_METHOD(data_long, "data-long", csv_table)
MAKE_READ_METHOD(group_rates, "group rates", csv_table)
MAKE_READ_METHOD(augm, "augmented", csv_table)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : print infer state  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize the output file with inferred parameters with column names and initial values.
void MCMC::print_infer_state(const std::string &output_filename) const
{
    if (hsp::na::is(output_filename))
        return;

    const Impl::Setup &setup = pImpl_->setup;
    const Impl::Parameter &param = pImpl_->param;
    const Impl::Augmented &augm = pImpl_->augm;

    std::ofstream ofs;
    ofs.open(output_filename.c_str());

    if (!ofs.is_open())
        throw std::runtime_error(
            "Unable to open file: " + output_filename);

    // Headers
    // header log-likelihood
    ofs << "iteration" << hsp::read::DELIMITER << "log_lik";

    // header parameters
    for (const std::string &p : param.list_sample)
        ofs << hsp::read::DELIMITER << p;

    // header augmented variable of id trace
    if (augm.is_any_draw())
        for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
            for (const std::string &a : augm.list_draws)
                ofs << hsp::read::DELIMITER << augm.smpl_individuals[k]->id() + "_" + a;

    // header chain
    ofs << hsp::read::DELIMITER << hsp::read::CHAIN_COLUMN << std::endl;

    // Values
    for (size_t c = 0; c < setup.n_chains; ++c)
        for (size_t k = 0; k < setup.c_sample_size; k++)
        {
            size_t i = k + c * setup.c_sample_size;

            // print iteration
            unsigned long sweep = setup.n_burnin + (k + 1) * setup.n_thin;
            ofs << sweep;

            // print log-likelihood
            ofs << hsp::read::DELIMITER << std::setprecision(2)
                << std::fixed << pImpl_->smpl_total_log_lik[i];

            // print parameter values
            for (const std::string &p : param.list_sample)
                ofs << hsp::read::DELIMITER << std::setprecision(5)
                    << std::fixed << param.smpl_values.find(p)->second[i];

            // print augmented variable of id trace
            if (!augm.smpl_values.empty())
                for (size_t k = 0; k < augm.smpl_individuals.size(); ++k)
                    for (const std::string &a : augm.list_draws)
                        ofs << hsp::read::DELIMITER << std::setprecision(5)
                            << std::fixed << augm.smpl_values[k].find(a)->second[i];

            // print which chain
            ofs << hsp::read::DELIMITER << (c + 1) << std::endl;
        }

    ofs.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC::Impl : print simulation  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Print all augmented variables of each individual.
std::ofstream open_ofstream_datawide(const std::string &f,
                                     const str_vec &data_headers,
                                     const str_vec &augm_headers)
{
    std::ofstream ofs;
    ofs.open(f.c_str());

    if (!ofs.is_open())
        throw std::runtime_error("Unable to open file: '" + f + "'");

    // print default data headers
    ofs << "id" << hsp::read::DELIMITER << "house_id" << hsp::read::DELIMITER << "status";

    // print other data headers
    for (const std::string &d : data_headers)
        ofs << hsp::read::DELIMITER << d;

    // print augmented variables headers
    for (const std::string &a : augm_headers)
        ofs << hsp::read::DELIMITER << a;

    ofs << std::endl;

    return ofs;
}

/// Print all augmented variables of each individual.
void update_ofstream_datawide(std::ofstream &ofs,
                              Household &house,
                              const str_vec &data_headers,
                              const str_vec &augm_headers)
{
    if (!ofs.is_open())
        throw std::runtime_error("Output file closed unexpectedly.");

    for (const Individual &member : house)
    {
        ofs << member.id()
            << hsp::read::DELIMITER << house.house_id()
            << hsp::read::DELIMITER << member.status();
        for (const std::string &d : data_headers)
            ofs << hsp::read::DELIMITER << member.data(d);

        for (const std::string &a : augm_headers)
            if (member.is_infected_status())
                ofs << hsp::read::DELIMITER << member.augm(a);
            else
                ofs << hsp::read::DELIMITER << hsp::na::DOUBLE;

        ofs << std::endl;
    }
}

/// Print the simulated data to the provided file.
void MCMC::print_simulated_data_wide(const std::string &pathtofile)
{
    Impl::Setup &setup = pImpl_->setup;
    const Impl::DataWide &data_wide = pImpl_->data_wide;
    const Impl::Augmented &augm = pImpl_->augm;

    try
    {
        // initialization
        setup.prog_declare("simulate");
        std::ofstream ofs = open_ofstream_datawide(pathtofile,
                                                   data_wide.headers,
                                                   augm.list_export);
        pImpl_->helper_init_simulation();
        // open progress bar
        setup.prog_open(data_wide.report_house.total);

        for (Household house : data_wide) // NOTE: deep copy of household
            try
            {
                // inplace simulation
                pImpl_->simulate_single_household_in_place(house);
                // print the simulated data
                update_ofstream_datawide(ofs, house,
                                         data_wide.headers,
                                         augm.list_export);
                // update progress bar
                setup.prog_update();
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While simulating " + house.err_name() + ":\n" + e.what());
            }

        // close progress bar
        setup.prog_close();

        // close ofstream
        if (!ofs.is_open())
            throw std::runtime_error(
                "Output file '" + pathtofile + "' closed unexpectedly.");
        ofs.close();
    }
    catch (const std::exception &e)
    {
        setup.prog_crash();
        hsp::msg::error(e.what());
    }
}

/// Estimate log_evidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
double MCMC::log_evidence(const std::string &scheme, unsigned long n_rungs)
{
    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = pImpl_->power_posteriors(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < pImpl_->setup.n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While estimating log_evidence with power posteriors:\n" +
                        std::string(e.what()));
    }
}
#endif
