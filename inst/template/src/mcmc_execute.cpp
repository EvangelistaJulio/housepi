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

#include <atomic>

#include "tools/string.h"
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
        setup.n_chains = 0;            // use 0 chains as a state variable
        assert_non_empty_data_wide_(); // check if data is empty
        warn_if_all_missing_status_(); // warn if all missing status
        param.user_reset();            // initialize non constant parameters
        augm.init();                   // initialize individual augmented relations
        augm.user_reset();             // initialize non constant augmented
        init_model_for_inference_();   // initialize model for inference
        update_full_loglike_(true);    // initialize cache + full loglikelihood
        return full_loglike_;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error(e.what());
    }

    return hsp::na::DOUBLE;
}
double MCMC::evaluate() { return pImpl_->evaluate(); }

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
    for (size_t R = 0; R < multiplicity.indexes.size(); ++R)
    {
        if (has_error.load(std::memory_order_relaxed))
            continue;

        Household &this_house = data_wide.restricted(R);
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
void MCMC::infer() { pImpl_->infer(); }

/// Assert setup definitions are ok to proceed with mcmc sampling.
void MCMC::Impl::init_monte_carlo_()
{
    setup.assert_mcmc_setup();     // assert if setup is ok for mcmc
    assert_non_empty_data_wide_(); // check if data is empty
    warn_if_all_missing_status_(); // warn if all missing status

    // check if there is actually something to sample
    if (!(param.is_any_sample() || augm.is_any_sample()))
        throw std::length_error(
            "Nothing to sample. Check your parameter and augmented definitions.");

    // initialize parameter lists for adaptation
    // NOTE: after all definitions finished to make cached parameters
    param.init_adaptable_lists();
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
        setup.prog_set_text();      // update progress bar text
        param.user_reset();         // reset parameter user definitions
        augm.user_reset();          // reset augmented user definitions
        update_full_loglike_(true); // init cache + full loglikelihood
        init_move_counts_();        // init move counts for adaptive

        // Finite adaptative iterations
        for (unsigned long i = 0; i < setup.n_adapt; i += setup.n_thin, setup.prog_update())
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
            {
                move_all_();
                adapt_rules_of_sampled_params_(i + t);
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
            update_posterior_sample_();
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

/// Adpat the rules of model parameters.
void MCMC::Impl::adapt_rules_of_sampled_params_(unsigned long step)
{
    const double ratio = static_cast<double>(step) / setup.n_adapt;

    // adpat all standard model parameters
    for (const std::string &p : param.list_adaptable_sample_standard)
        param.at(p).adapt_rule(ratio);

    // adapt all parameters with cached dependecies
    for (const std::string &p : param.list_adaptable_sample_cached)
        param.at(p).adapt_rule(ratio);
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
    for (Household &house : data_wide.restricted())
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
    smpl_values_.resize(pImpl_.multiplicity.indexes.size());
    for (size_t R = 0; R < pImpl_.multiplicity.indexes.size(); ++R)
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
    for (size_t R = 0; R < pImpl_.multiplicity.indexes.size(); ++R)
    {
        Household &house = house_vector_[pImpl_.multiplicity.indexes[R]];
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
        data_wide.init_draws(sample_size);
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
        data_wide.update_draws(smpl_count_);
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

    data_wide.clear_draws();
    param.clear_draws();
    augm.clear_draws();
}
void MCMC::clear_memory() { pImpl_->clear_memory(); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : estimate logevidence  /////////////////////////////
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
std::vector<dbl_vec> MCMC::Impl::power_posteriors_(const std::string &scheme, unsigned long n_rungs)
{
    // check scheme choice
    const str_vec options = {"adaptive", "power-fractions"};
    if (!(options % scheme))
        throw std::invalid_argument(
            "Invoked invalid power posterior scheme \"" + scheme + "\"." +
            hsp::str::find_closest(scheme, options));

    // check number of rungs for the thermodynamic ladder
    if (n_rungs < 2)
        throw std::invalid_argument(
            "N. of rungs was not greater than 1.");
    setup.n_rungs = n_rungs; // NOTE: set before progress bar

    // initialize mcmc sampling
    init_monte_carlo_();

    // buffer user setup
    double buffer_inv_T = setup.inv_T;

    setup.prog_declare();        // declare current progress bar
    init_model_for_inference_(); // initialize model for inference
    augm.init();                 // initialize augmented variables
    setup.prog_open((setup.c_n_sweeps * setup.n_chains) * (setup.n_rungs + 1));

    // thermodynamic ladder estimates
    dbl_vec itm = {}; // inverse temperature
    dbl_vec avg = {}; // mean loglikelihood
    dbl_vec var = {}; // variance loglikelihood

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
    run_monte_carlo_();
    std::pair<double, double> stats = compute_stats(smpl_full_loglike_);

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
    run_monte_carlo_();
    std::pair<double, double> stats = compute_stats(smpl_full_loglike_);

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
        if ((new_t <= itm[k]) || (new_t >= itm[k + 1]))
            new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);

        // power posterior sampling
        setup.inv_T = new_t;
        run_monte_carlo_();
        stats = compute_stats(smpl_full_loglike_);

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
        run_monte_carlo_();

        std::pair<double, double> stats = compute_stats(smpl_full_loglike_);
        itm.push_back(setup.inv_T);
        avg.push_back(stats.first);
        var.push_back(stats.second);
    }
}

#ifdef _HOUSEPI_RCPP
/// Estimate logevidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
Rcpp::List MCMC::Impl::logevidence(const std::string &scheme, unsigned long n_rungs)
{
    using namespace Rcpp;

    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = power_posteriors_(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < setup.n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return List::create(_["estimate"] = output,
                            _["scheme"] = scheme,
                            _["n_rungs"] = setup.n_rungs,
                            _["inv_temper"] = wrap(inv_tem),
                            _["mean_loglike"] = wrap(avg_vec),
                            _["var_loglike"] = wrap(var_vec));
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While estimating logevidence with power posteriors:\n" +
                        std::string(e.what()));
    }
    return List();
}
Rcpp::List MCMC::logevidence(const std::string &scheme, unsigned long n_rungs) { return pImpl_->logevidence(scheme, n_rungs); }

#else
/// Estimate logevidence according to a given scheme and number of rungs.
/// Uses trapezoidal rule to estimate logarithm of marginal likelihood (evidence).
///
/// Ref.
/// -----
/// Ref: Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
double MCMC::Impl::logevidence(const std::string &scheme, unsigned long n_rungs)
{
    try
    {
        double output = 0;
        const std::vector<dbl_vec> pp = power_posteriors_(scheme, n_rungs);
        const dbl_vec &inv_tem = pp[0];
        const dbl_vec &avg_vec = pp[1];
        const dbl_vec &var_vec = pp[2];

        for (unsigned long k = 0; k < setup.n_rungs; ++k)
        {
            double d_inv_tem = inv_tem[k + 1] - inv_tem[k];
            output += (avg_vec[k + 1] + avg_vec[k]) * d_inv_tem / 2.;
            output += (var_vec[k] - var_vec[k + 1]) * d_inv_tem * d_inv_tem / 12.;
        }

        return output;
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While estimating logevidence with power posteriors:\n" +
                        std::string(e.what()));
    }

    return hsp::na::DOUBLE;
}
double MCMC::logevidence(const std::string &scheme, unsigned long n_rungs) { return pImpl_->logevidence(scheme, n_rungs); }
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMC : Impl : simulation  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _HOUSEPI_RCPP
/// Perform a unique repetition.
Rcpp::DataFrame MCMC::Impl::simulate_repetition()
{
    using namespace Rcpp;

    const size_t n_rows = data_wide.report_indiv.total;

    IntegerVector df_status = IntegerVector(n_rows);

    std::vector<NumericVector> df_augm;
    df_augm.reserve(augm.list_export.size());
    for (size_t c = 0; c < augm.list_export.size(); ++c)
        df_augm.emplace_back(n_rows);

    size_t r = 0;
    for (Household house : data_wide) // NOTE: deep copy of household
        try
        {
            simulate_single_household_(house); // inplace simulation
            for (const Individual &member : house)
            {
                df_status[r] = member.status();
                for (size_t c = 0; c < augm.list_export.size(); ++c)
                    df_augm[c][r] = hsp_to_rcpp(member.augm(augm.list_export[c]));
                ++r;
            }
        }
        catch (const std::exception &e)
        {
            throw std::runtime_error(
                "While simulating " + house.err_name() + ":\n" +
                std::string(e.what()));
        }

    DataFrame output = DataFrame::create(_["status"] = df_status);
    for (size_t c = 0; c < augm.list_export.size(); ++c)
        output.push_back(df_augm[c], augm.list_export[c]);

    return output;
}

/// Perform simulation and return simulated data.
Rcpp::List MCMC::Impl::simulate(size_t n_rep)
{
    using namespace Rcpp;
    std::atomic<bool> has_error(false);

    try
    {
        // initialization
        aux_init_simulate_();
        setup.prog_open(data_wide.report_house.total * n_rep);

        // run simulations
        std::vector<DataFrame> results(n_rep);

#pragma omp parallel for
        for (size_t rep = 0; rep < n_rep; ++rep)
        {
            if (has_error.load(std::memory_order_relaxed))
                continue;
            try
            {
                results[rep] = simulate_repetition();
            }
            catch (const std::exception &e)
            {
#pragma omp critical
                {
                    if (!has_error.exchange(true))
                        error_message = e.what();
                }
            }
        }
        if (!error_message.empty())
            throw std::runtime_error(error_message);

        // convert to Rcpp List
        List output(n_rep);
        for (size_t rep = 0; rep < n_rep; ++rep)
            output[rep] = std::move(results[rep]);

        // close progress bar
        setup.prog_close();

        return output;
    }
    catch (std::exception &e)
    {
        hsp::msg::error(e.what());
    }

    return List::create();
}
Rcpp::List MCMC::simulate(size_t n_rep) { return pImpl_->simulate(n_rep); }
#endif
