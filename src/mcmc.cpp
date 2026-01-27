/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <algorithm>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "tools/atomic.h"
#include "tools/keys_helper.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/string.h"
#include "tools/types.h"

#include "augmented.h"
#include "data.h"
#include "data_reports.h"
#include "group_rate.h"
#include "infec_profile.h"
#include "mcmc_impl.h"
#include "mcmc.h" /// UPDATE: should remove I/O and Rcpp dependencies from MCMC...
#include "mcmc_impl.h"
#include "parameter.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  *** MCMCImpl : clas ***  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

MCMCImpl::MCMCImpl(
    const std::vector<std::string> &data_wide_keys,
    const std::vector<std::string> &hh_size_t_keys,
    const std::vector<std::string> &hh_double_keys,
    const std::vector<std::string> &hh_matrix_keys,
    const std::vector<std::string> &data_long_keys,
    const std::vector<std::string> &param_keys,
    const std::vector<std::string> &group_rate_keys,
    const std::vector<std::string> &augm_keys,
    std::function<void()> func_init_inference,
    std::function<double(Household &)> func_household_log_lik,
    std::function<void(Household &, Individual &)> func_update_cache)
    : atomic_error(std::make_unique<hsp::atomic::Error>()),
      // modules
      setup(std::make_shared<Setup>()),
      data_wide(std::make_shared<DataWide>(data_wide_keys, hh_size_t_keys, hh_double_keys, hh_matrix_keys)),
      data_long(std::make_shared<DataLong>(data_long_keys, data_wide)),
      param(std::make_shared<Parameter>(param_keys)),
      group_rate(std::make_shared<GroupRate>(group_rate_keys, data_wide, param)),
      infec_profile(std::make_shared<InfecProfile>(param)),
      augm(std::make_shared<Augmented>(augm_keys, data_wide, param, infec_profile)),
      // callbacks
      callback_init_inference_(std::move(func_init_inference)),
      callback_house_log_lik_(std::move(func_household_log_lik)),
      callback_update_cache_(std::move(func_update_cache))
{
    /// UPDATE: waiting for refactoring of read_attribute
    param->link_to_data_wide(data_wide->pRegistry);
    param->link_to_group_rate(group_rate->pRegistry);
    param->link_to_augm(augm->pRegistry);
}
MCMCImpl::~MCMCImpl() = default;

void MCMCImpl::init_household_cache()
{
    // init cache dependencies for all households and individuals
    for (size_t R = 0; R < data_wide->head_class_house_indexes.size(); ++R)
        data_wide->init_house_cache(data_wide->restricted(R));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMCImpl : evaluate log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assert setup definitions are ok to proceed with log-likelihood evaluation.
void MCMCImpl::init_evaluate()
{
    setup->n_chains = 0;                  // use 0 chains as a state variable
    data_wide->assert_non_empty();        // check if data is empty
    data_wide->warn_all_missing_status(); // warn if all missing status
    callback_init_inference_();           // user initialization of model for inference
    init_household_cache();               // cache user household variables
    augm->init();                         // initialize individual augmented relations
    param->user_reset(setup->prng());     // initialize non constant parameters
    augm->user_reset(setup->prng());      // initialize non constant augmented
}

/// Evaluate log-likelihood for the given household.
double MCMCImpl::evaluate_single_house_log_lik(std::string_view house_id)
{
    atomic_error->restore(); // reset atomic error
    infec_profile->update(); // eventually buffer infectivity profile
    double output = atomic_materialize_house_log_lik_(data_wide->at(house_id));
    atomic_error->throw_exception(); // catch and throw if previous failed
    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMCImpl : update house log-likelihood  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute new log-likelihood, store old value in buffer. Same for households.
double MCMCImpl::update_total_log_lik(bool is_cache_depen)
{
    double new_ll = 0.0;

    // buffer infectivity profile if lookup tables are active
    if (is_cache_depen)
        infec_profile->update();

    atomic_error->restore();
#pragma omp parallel for reduction(+ : new_ll)
    for (size_t R = 0; R < data_wide->head_class_house_indexes.size(); ++R)
    {
        if (atomic_error->is_pending())
            continue;

        new_ll += atomic_materialize_house_log_lik_(data_wide->restricted(R), is_cache_depen);
    }
    atomic_error->throw_exception();

    return new_ll;
}

/// Update and return contribution of given household to full logliklehood.
/// NOTE: house multiplicity doesn't work with data augmentation
double MCMCImpl::atomic_materialize_house_log_lik_(Household &house, bool is_cache_depen)
{
    double output = 0;
    try
    {
        // buffer all household cache dependencies
        if (is_cache_depen)
        {
            house.buffer_all_cache();
            for (auto &member : house.infected())
                callback_update_cache_(house, member);
        }

        // compute and set household observation log-likelihood
        output = callback_house_log_lik_(house);
        house.set_obs_log_lik(output);

        // compute and set individual augmented log-likelihood
        if (augm->is_any_model())
            for (auto &member : house.infected())
                output += augm->materialize_augm_model_log_lik(member);
    }
    catch (const std::exception &e)
    {
        atomic_error->set_if_first("While computing log-likelihood for " +
                                   house.err_name() + ":\n" + std::string(e.what()));
    }

    return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMCImpl : inference  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assert setup definitions are ok to proceed with mcmc sampling.
void MCMCImpl::init_monte_carlo()
{
    // check if data is empty
    data_wide->assert_non_empty();

    // check if there is actually something to sample
    if (!param->is_any_sample() && !augm->is_any_sample())
        throw std::length_error(
            "Nothing to sample. Check your parameter and augmented definitions.");
    // warn if all missing status
    data_wide->warn_all_missing_status();

    // assert if setup is ok for mcmc
    setup->assert_mcmc_setup();

    // Drop adaptation steps if nothing to adapt
    // UPDATE: if adaptation applied to augmented variables, update this
    if (param->list_adaptable_sample_cached.empty() &&
        param->list_adaptable_sample_standard.empty())
        setup->drop_mcmc_adapt();

    // Initialize individual augmented sample
    augm->init_individual_sample(setup->id_draws);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMCImpl : run monte carlo  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Markov chain core, used for inference and power posteriors.
/// NOTE: sweep is the number of previously completed thin sweeps.
void MCMCImpl::run_monte_carlo()
{
    init_posterior_sample_(); // initialize posterior sample

    // run multiple MCMC chains&
    for (setup->chain = 1; setup->chain <= setup->n_chains; ++setup->chain)
    {
        init_monte_carlo_chain_();

        // Finite adaptative iterations
        for (unsigned long i = 0; i < setup->n_adapt; i += setup->n_thin, setup->prog_update())
            for (unsigned long t = 1; t <= setup->n_thin; ++t)
            {
                move_all_();
                adapt_proposal_of_sampled_params_(i + t);
            }

        // Burn-in iterations
        for (unsigned long i = 0; i < setup->n_burnin; i += setup->n_thin, setup->prog_update())
            for (unsigned long t = 1; t <= setup->n_thin; ++t)
                move_all_();

        // Posterior sampling
        for (unsigned long i = setup->n_burnin; i < setup->n_steps; i += setup->n_thin, setup->prog_update())
        {
            for (unsigned long t = 1; t <= setup->n_thin; ++t)
                move_all_();
            store_posterior_sample_();
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMCImpl : monte carlo moves /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Move all model parameters and augmented variables.
void MCMCImpl::move_all_()
{
    // move augmented variables to sample, if any
    if (augm->is_any_sample())
        move_all_augm_();

    // move all standard model parameters
    for (const auto &p : param->list_adaptable_sample_standard)
        param->at(p).update_accept_rate(move_single_parameter_(p, false));
    for (const auto &p : param->list_non_adapt_sample_standard)
        move_single_parameter_(p, false);

    // move parameters with cached dependecies
    for (const auto &p : param->list_adaptable_sample_cached)
        param->at(p).update_accept_rate(move_single_parameter_(p, true));
    for (const auto &p : param->list_non_adapt_sample_cached)
        move_single_parameter_(p, true);
}

/// Attempt to move a single parameter for the given inverse temperature.
bool MCMCImpl::move_single_parameter_(size_t p, bool is_cache_depen)
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
void MCMCImpl::adapt_proposal_of_sampled_params_(unsigned long step)
{
    const double ratio = static_cast<double>(step) / static_cast<double>(setup->n_adapt);

    // adpat all standard model parameters
    for (const auto &p : param->list_adaptable_sample_standard)
        param->at(p).adapt_proposal(ratio);

    // adapt all parameters with cached dependecies
    for (const auto &p : param->list_adaptable_sample_cached)
        param->at(p).adapt_proposal(ratio);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMCImpl : parameter moves  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Intialize a single Markov chain for Monte Carlo.
void MCMCImpl::init_monte_carlo_chain_()
{
    setup->prog_set_text();                     // set progress bar text
    param->user_reset(setup->prng());           // reset parameter user definitions
    augm->user_reset(setup->prng());            // reset augmented user definitions
    set_total_log_lik_(update_total_log_lik()); // init cached, set total log-likelihood

    // init move counts for adaptive
    for (const auto p : param->list_sample)
        param->at(p).reset_accept_rate();
}

/// propose a valid move (return true), otherwise update refusal (return false).
bool MCMCImpl::propose_valid_param_move_(size_t p)
{
    auto &this_par = param->at(p);

    // set new value of the parameter and store the old to buffer
    this_par.propose(setup->prng());

    // recover buffer if not valid move
    if (!this_par.is_valid_move())
    {
        this_par.recover_value();
        return false;
    }

    return true;
}

/// Restore the parameter value in the buffer.
void MCMCImpl::recover_parameter_(size_t p, bool is_cache_depen)
{
    // recover buffer value of parameter
    param->at(p).recover_value();

    // recover buffer value of full log-like
    total_log_lik_ = total_buffer_;

    // recover buffered log-likelihood in households
    for (size_t R = 0; R < data_wide->head_class_house_indexes.size(); ++R)
    {
        auto &house = data_wide->restricted(R);
        house.recover_buffer_values(is_cache_depen);
        if (augm->is_any_sample())
            for (auto &member : house.infected())
                member.recover_augm_model_log_lik();
    }

    if (is_cache_depen)
        infec_profile->recover();
}

/// Metropolis-Hastings to accept or not current parameter move.
bool MCMCImpl::accept_valid_param_move_(size_t p)
{
    const double D_ll = total_log_lik_ - total_buffer_;
    const double move_logratio = param->at(p).get_move_logratio();
    return setup->unif() < std::exp(setup->inv_T * D_ll + move_logratio);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  MCMCImpl : move all augmented variables  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Move augmented variables from all households.
void MCMCImpl::move_all_augm_()
{
    double delta_ll = 0.;
    atomic_error->restore();

#pragma omp parallel for reduction(+ : delta_ll)
    for (size_t H = 0; H < data_wide->report_house->total; ++H)
    {
        if (atomic_error->is_pending())
            continue;

        delta_ll += atomic_move_augm_single_house_(H);
    }
    atomic_error->throw_exception();

    // increment total log-likelohood by a delta, buffer old value
    set_total_log_lik_(total_log_lik_ + delta_ll);
}

/// Move augmented variables of infectious members inside household.
/// The move is performed one member at a time, one augmented variable at a time.
double MCMCImpl::atomic_move_augm_single_house_(size_t H)
{
    auto &house = data_wide->at(H);
    double delta_ll = 0.;

    try
    {
        // loop over infectious residents of the househould
        for (auto &member_i : house.infected())
            // loop over augmented variables of infectious
            for (auto a : augm->status_to_sample[member_i.status()])
                delta_ll += metropolis_augm_single_move_(house, member_i, a);
    }
    catch (const std::exception &e)
    {
        atomic_error->set_if_first("While sampling augmented variables for " +
                                   house.err_name() + ":\n" + std::string(e.what()));
    }

    return delta_ll;
}

/// Take care of a a single augmented Metropolis move.
double MCMCImpl::metropolis_augm_single_move_(Household &house, Individual &member, size_t a)
{
    auto &xaugm = member.get_augm_variable(a);
    xaugm.propose(setup->prng()); // propose new value and buffer old

    // proceed only if the new value is consistent
    if (xaugm.is_valid_move())
    {
        // update dependencies to perform computations
        house.buffer_all_cache();
        callback_update_cache_(house, member);

        // Metropolis-Hastings
        // NOTE: materialize calls must come after getters
        double delta_ll = -house.get_obs_log_lik();
        delta_ll += augm->materialize_augm_model_logratio(member, a);
        // compute and set household observation log-likelihood
        const double house_ll = callback_house_log_lik_(house);
        house.set_obs_log_lik(house_ll);

        delta_ll += house_ll;

        if (setup->unif() < std::exp(setup->inv_T * delta_ll - xaugm.get_propose_logratio()))
            return delta_ll;

        // else recover cached values
        house.recover_buffer_values(true);
        member.recover_augm_model_log_lik();
    }
    // else recover buffered value
    xaugm.recover_value();

    return 0.;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMCImpl : posterior sample  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Init posterior sampling.
void MCMCImpl::init_posterior_sample_()
{
    unsigned long sample_size = setup->c_sample_size * setup->n_chains;

    // diagnostics over full log-likelihood
    smpl_count = 0;
    smpl_total_log_lik.clear();
    smpl_total_log_lik.reserve(sample_size);
    smpl_total_log_lik.assign(sample_size, hsp::na::DOUBLE);

    // initialize  other diagnostic vectors
    // UPDATE: I need to store also power posterior sample
    // NOTE: it will depend on chosen rungs in that case
    if (!setup->is_power_post())
    {
        data_wide->init_draws(sample_size);
        param->init_draws(sample_size);
        augm->init_draws(sample_size);
    }
}

/// Store current Bayesian posterior sample.
void MCMCImpl::store_posterior_sample_()
{
    // store sample of the full log-likelihood
    smpl_total_log_lik[smpl_count] = total_log_lik_;

    if (!setup->is_power_post())
    {
        data_wide->store_draws(smpl_count);
        param->store_draws(smpl_count);
        augm->store_draws(smpl_count);
    }

    ++smpl_count;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMCImpl : estimate log_evidence  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace
{
    /// Welford online's algorithm for mean and variance.
    /// Ref :
    /// ----
    /// https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    std::pair<double, double> compute_stats(const std::vector<double> &data)
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

        return {avg, sqr_diff / static_cast<double>(N - 1.)};
    }
}

/// Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
///
/// Ref
/// -----
/// Friel, N., Hurn, M. & Wyse, J.
/// Improving power posterior estimation of statistical evidence.
/// Stat Comput 24, 709–723 (2014).
std::vector<std::vector<double>> MCMCImpl::power_posteriors(std::string_view scheme, unsigned long n_rungs)
{
    // check scheme choice
    const std::vector<std::string> options = {"adaptive", "power-fractions"};
    if (std::find(options.begin(), options.end(), scheme) == options.end())
        throw std::invalid_argument(
            "Invoked invalid power posterior scheme '" + std::string(scheme) + "'." +
            hsp::str::find_closest(scheme, options));

    // check number of rungs for the thermodynamic ladder
    if (n_rungs < 2)
        throw std::invalid_argument("N. of rungs was not greater than 1.");

    setup->n_rungs = n_rungs;           // NOTE: set before progress bar
    double buffer_inv_T = setup->inv_T; // buffer user setup
    init_monte_carlo();                 // initialize mcmc sampling
    callback_init_inference_();         // user initialization of model for inference
    init_household_cache();             // cache user household variables
    augm->init();                       // initialize augmented variables
    setup->prog_declare();              // declare current progress bar
    setup->prog_open_for_power_post();  // print progress bar

    // thermodynamic ladder estimates
    std::vector<double> itm = {}; // inverse temperature
    std::vector<double> avg = {}; // mean log-likelihood
    std::vector<double> var = {}; // variance log-likelihood

    // run power posteriors
    power_post_zero_(itm, avg, var);
    if (scheme == "adaptive")
        power_post_adaptive_(itm, avg, var);
    else if (scheme == "power-fractions")
        power_post_fractions_(itm, avg, var);

    setup->prog_close();         // close progress bar
    setup->inv_T = buffer_inv_T; // restore user setup

    return {itm, avg, var};
}

/// Power posteriors at 0 inverse temperature.
void MCMCImpl::power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
{
    // UPDATE: should I use previous inference for itm = 1 if available?
    // UPDATE: always under-estimate time if starting from itm = 0

    setup->rung = 0;
    setup->inv_T = 0.0;
    run_monte_carlo();
    const auto [stats_avg, stats_var] = compute_stats(smpl_total_log_lik);

    itm.push_back(setup->inv_T);
    avg.push_back(stats_avg);
    var.push_back(stats_var);
}

/// Power posteriors with adaptive scheme for the thermodynamic ladder.
void MCMCImpl::power_post_adaptive_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
{
    // Inverse Temperature 1
    {
        setup->rung = 1;
        setup->inv_T = 1.0;
        run_monte_carlo();
        const auto [stats_avg, stats_var] = compute_stats(smpl_total_log_lik);

        itm.push_back(setup->inv_T);
        avg.push_back(stats_avg);
        var.push_back(stats_var);
    }

    // Inverse Temperatures (0, 1)
    size_t k = 0; // index to interval lower bound
    bool is_right = false;
    while (setup->rung < setup->n_rungs)
    {
        ++setup->rung;

        // adaptive choice of inverse temperature
        // intersection between tangents
        double new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
        // otherwise weighted average
        if (new_t <= itm[k] || new_t >= itm[k + 1])
            new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);

        // power posterior sampling
        setup->inv_T = new_t;
        run_monte_carlo();
        const auto [stats_avg, stats_var] = compute_stats(smpl_total_log_lik);

        // find larger rectangular discretization error
        is_right = (itm[k + 1] - setup->inv_T) * (avg[k + 1] - stats_avg) > (new_t - itm[k]) * (stats_avg - avg[k]);

        // insert in-between for all lists
        itm.insert(itm.begin() + k + 1, setup->inv_T);
        avg.insert(avg.begin() + k + 1, stats_avg);
        var.insert(var.begin() + k + 1, stats_var);

        // if right interval is chosen, index points to inserted as lower bound
        k += (is_right);
    }
}

/// Power posteriors with power fractions scheme for the thermodynamic ladder.
void MCMCImpl::power_post_fractions_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
{
    for (setup->rung = 1; setup->rung <= setup->n_rungs; ++setup->rung)
    {
        // inverse tempearture = (rung / #rungs)^5
        setup->inv_T = std::pow(static_cast<double>(setup->rung) / static_cast<double>(setup->n_rungs), 5.);
        run_monte_carlo();

        const auto [stats_avg, stats_var] = compute_stats(smpl_total_log_lik);
        itm.push_back(setup->inv_T);
        avg.push_back(stats_avg);
        var.push_back(stats_var);
    }
}

/// Auxiliary method to initialize simulations.
/// UPDATE: once I allow for multi-treading and a flag other than n_chains = 0, I can drop this.
void MCMCImpl::helper_init_simulation()
{
    data_wide->assert_non_empty(); // check if data is empty
    setup->n_chains = 0;           // use 0 chains as a state variable
    setup->n_threads = 1;          // impose no parallelization
    setup->set_n_threads();        // change accordingly
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : constructor  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::MCMC(
    const std::vector<std::string> &data_wide_keys,
    const std::vector<std::string> &hh_size_t_keys,
    const std::vector<std::string> &hh_double_keys,
    const std::vector<std::string> &hh_matrix_keys,
    const std::vector<std::string> &data_long_keys,
    const std::vector<std::string> &param_keys,
    const std::vector<std::string> &group_rate_keys,
    const std::vector<std::string> &augm_keys)
    : pImpl_(std::make_unique<MCMCImpl>(
          data_wide_keys, hh_size_t_keys, hh_double_keys, hh_matrix_keys, data_long_keys,
          param_keys, group_rate_keys, augm_keys,
          [this]()
          { this->init_model_for_inference(); },
          [this](Household &h) -> double
          { return this->compute_house_obs_log_lik(h); },
          [this](Household &h, Individual &i)
          { this->update_cache_dependecies(h, i); })) {}
MCMC::~MCMC() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMC : general methods  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double MCMC::uniform() { return pImpl_->setup->unif(); }
double MCMC::par(size_t p) const { return pImpl_->param->at(p).get(); }

size_t MCMC::get_max_hh_cases() const { return pImpl_->data_wide->report_hh_size->max_cases.first; }
void MCMC::assert_no_missing_status() const
{
    if (pImpl_->data_wide->report_indiv->missing > 0)
        hsp::msg::error(
            "Individual with missing status are incompatible with current model.");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : contact rates ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Build all group rates and contact rates and set rate per group. // UPDATE
void MCMC::build_all_rates()
{
    try
    {
        for (auto &house : *pImpl_->data_wide)
            init_contact_rates_(house);
    }
    catch (const std::exception &e)
    {
        hsp::msg::error("While initializing contact rates:\n" +
                        std::string(e.what()));
    }
}

/// Set contact rate according to user assignement. // UPDATE
void MCMC::init_contact_rates_(Household &house)
{
    // loop over all couples
    for (const auto &member_k : house)     // from : member_k
        for (const auto &member_i : house) // -> to : member_i
            if (member_k != member_i)
            {
                pair_contact_rate_ = "";
                init_house_contact_rates(house, member_k, member_i);

                if (!hsp::na::is(pair_contact_rate_))
                {
                    // default reference rate
                    Attribute y(hsp::rates::VALUE_REFERENCE);

                    if (pair_contact_rate_ == hsp::rates::MISSING)
                        throw std::invalid_argument(
                            "Missing contact rate for individuals '" +
                            member_k.id() + "' and '" + member_i.id() +
                            "' (" + house.err_name() + ").");

                    else if (pair_contact_rate_ != hsp::rates::REFERENCE) // parameter rate
                        y = pImpl_->param->at(pImpl_->param->pRegistry->from_string(pair_contact_rate_));

                    house.init_contact_rates(static_cast<size_t>(member_k), static_cast<size_t>(member_i), y);
                }
            }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  MCMC : cache attributes  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void MCMC::init_infec_profile_template(std::vector<std::pair<std::string, std::string>> x)
{
    hsp::str_table input = {};
    for (const auto &[key, val] : x)
        input[key] = {val};
    pImpl_->infec_profile->init_template(std::move(input));
}
/// Request discrete infectivity profile.
void MCMC::request_infec_profile_discrete() { pImpl_->infec_profile->request_discrete(); }
double MCMC::infec_profile_inst(const Individual &i, double t) const
{
    return pImpl_->infec_profile->inst(&i, t);
}
double MCMC::infec_profile_cumul(const Individual &i, double t) const
{
    return pImpl_->infec_profile->cumul(&i, t);
}

void MCMC::clear_all_augm(Individual &member) { pImpl_->augm->clear_all(member); }
void MCMC::init_missing_augm(Individual &member)
{
    pImpl_->augm->init_missing(member);
    pImpl_->augm->user_reset(pImpl_->setup->prng(), member);
}

void MCMC::request_discrete_augm(size_t a) { pImpl_->augm->user_discrete_set.emplace(a); }
void MCMC::request_discrete_param(size_t p) { pImpl_->param->user_discrete_set.emplace(p); }

/// Initialize model and evaluate log-likelihood.
double MCMC::evaluate_log_lik(std::string_view house_id)
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

/// Performs the inference using Metropolis-Hastings.
void MCMC::run_inference()
{
    try
    {
        pImpl_->init_monte_carlo();               // initialize mcmc sampling
        init_model_for_inference();               // user initialization model
        pImpl_->init_household_cache();           // cache user household variables
        pImpl_->augm->init();                     // initialize augmented variables
        pImpl_->setup->prog_declare();            // declare current progress bar
        pImpl_->setup->prog_open_for_inference(); // open progress bar
        pImpl_->run_monte_carlo();                // run Monte Carlo
        pImpl_->setup->prog_close();              // close progress bar
    }
    catch (std::exception &e)
    {
        pImpl_->setup->prog_crash();
        hsp::msg::error(e.what());
    }
}
