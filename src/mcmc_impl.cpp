/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/atomic.h"
#include "tools/keys.h"
#include "tools/loglike.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/openmp.h"
#include "tools/string.h"
#include "tools/types.h"

#include "mcmc.h"
#include "data_wide.h"
#include "household.h"
#include "data_long.h"
#include "parameter.h"
#include "group_rate.h"
#include "contact_rate.h"
#include "infec_profile.h"
#include "augmented.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  *** MCMCBase : class ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

MCMCBase::MCMCBase()
    : data_long(&data_wide),
      group_rate(&data_wide, &param),
      infec_profile(&param),
      augm(&data_wide, &param, &infec_profile) {}
MCMCBase::~MCMCBase() = default;

void MCMCBase::init_registry(const std::vector<std::string> &data_wide_keys,
                             const std::vector<std::string> &hh_size_keys,
                             const std::vector<std::string> &hh_size_vec_keys,
                             const std::vector<std::string> &hh_double_vec_keys,
                             const std::vector<std::string> &hh_double_mtx_keys,
                             const std::vector<std::string> &data_long_keys,
                             const std::vector<std::string> &param_keys,
                             const std::vector<std::string> &group_rate_keys,
                             const std::vector<std::string> &augm_keys)
{
    data_wide.init_registry(data_wide_keys,
                            hh_size_keys,
                            hh_size_vec_keys,
                            hh_double_vec_keys,
                            hh_double_mtx_keys);
    data_long.init_registry(data_long_keys);
    param.init_registry(param_keys);
    group_rate.init_registry(group_rate_keys);
    augm.init_registry(augm_keys);

    /// UPDATE: waiting for refactoring of read_attribute
    param.bind_to_data_wide(data_wide.pRegistry.get());
    param.bind_to_group_rate(group_rate.pRegistry.get());
    param.bind_to_augm(augm.pRegistry.get());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  MCMCBase : monte carlo  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMCBase::init_monte_carlo_()
{
    // check if there is actually something to sample
    if (!param.is_any_sample() && !augm.is_any_sample())
        throw std::length_error(
            "Nothing to sample. Check your parameter and augmented definitions.");

    // assert if setup is ok for mcmc
    setup.assert_mcmc_setup();

    // Drop adaptation steps if nothing to adapt
    // UPDATE: if adaptation applied to augmented variables, update this
    if (param.list_adaptable_sample_cached.empty() &&
        param.list_adaptable_sample_standard.empty())
        setup.drop_mcmc_adapt();

    // Initialize individual augmented sample
    augm.init_individual_sample(setup.id_draws);
}

void MCMCBase::init_monte_carlo_chain_()
{
    setup.prog_set_text();          // set progress bar text
    param.user_reset(setup.prng()); // reset parameter user definitions
    augm.user_reset(setup.prng());  // reset augmented user definitions

    // init move counts for adaptive
    for (size_t p : param.list_sample)
        param.at(p).reset_accept_rate();
}

void MCMCBase::recover_parameter_(size_t p, bool is_cache_depen)
{
    // recover buffer value of parameter
    param.at(p).recover();
    // recover buffer value of full log-like
    tot_log_lik.recover();

    // recover buffered log-likelihood in households
    for (size_t R = 0; R < data_wide.restricted_size(); ++R)
    {
        Household &house = data_wide.restricted(R);
        house.obs_log_lik_.recover();

        if (augm.is_any_sample())
            for (size_t i : house.infected_hhindexes())
                house.member(i).augm_log_lik_.recover();
    }

    if (is_cache_depen)
    {
        data_wide.recover_all_buffer();
        infec_profile.recover();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  MCMCBase : posterior sample  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void MCMCBase::init_posterior_sample_()
{
    unsigned long sample_size = setup.c_sample_size * setup.n_chains;

    if (!setup.is_power_post())
    {
        data_wide.init_draws(sample_size);
        param.init_draws(sample_size);
        augm.init_draws(sample_size);
    }
    tot_log_lik.init_draws(sample_size);
}

void MCMCBase::store_posterior_sample_()
{
    if (!setup.is_power_post())
    {
        data_wide.store_draws(tot_log_lik.count());
        param.store_draws(tot_log_lik.count());
        augm.store_draws(tot_log_lik.count());
    }
    // store total log-likelihood and update count
    tot_log_lik.store_draws();
}

/* --FUTURE-LOG-EVIDENCE--
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  MCMC : estimate log_evidence  ////////////////////////////////
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

std::vector<std::vector<double>> MCMC::power_posteriors(std::string_view scheme, unsigned long n_rungs)
{
    // check scheme choice
    const std::vector<std::string> options = {"adaptive", "power-fractions"};
    if (std::find(options.begin(), options.end(), scheme) == options.end())
        throw std::invalid_argument(
            "Invoked invalid power posterior scheme '" + std::string(scheme) + "'." +
            hsp::str::find_closest(scheme, options));

    // check number of rungs for the thermodynamic ladder
    if (n_rungs < 2)
        throw std::invalid_argument(
            "N. of rungs was not greater than 1.");

    setup.n_rungs = n_rungs;           // NOTE: set before progress bar
    double buffer_inv_T = setup.inv_T; // buffer user setup
    init_evaluate_(false);             // initialize for evaluation
    init_monte_carlo_();               // initialize mcmc sampling
    setup.prog_declare();              // declare current progress bar
    setup.prog_open_for_power_post();  // print progress bar

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

    setup.prog_close();         // close progress bar
    setup.inv_T = buffer_inv_T; // restore user setup

    return {itm, avg, var};
}

void MCMC::power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
{
    // UPDATE: should I use previous inference for itm = 1 if available?
    // UPDATE: always under-estimate time if starting from itm = 0

    setup.rung = 0;
    setup.inv_T = 0.0;
    run_monte_carlo_();
    const auto [stats_avg, stats_var] = compute_stats(tot_log_lik.smpl_values);

    itm.push_back(setup.inv_T);
    avg.push_back(stats_avg);
    var.push_back(stats_var);
}

void MCMC::power_post_adaptive_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
{
    // Inverse Temperature 1
    {
        setup.rung = 1;
        setup.inv_T = 1.0;
        run_monte_carlo_();
        const auto [stats_avg, stats_var] = compute_stats(tot_log_lik.smpl_values);

        itm.push_back(setup.inv_T);
        avg.push_back(stats_avg);
        var.push_back(stats_var);
    }

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
        run_monte_carlo_();
        const auto [stats_avg, stats_var] = compute_stats(tot_log_lik.smpl_values);

        // find larger rectangular discretization error
        is_right = (itm[k + 1] - setup.inv_T) * (avg[k + 1] - stats_avg) > (new_t - itm[k]) * (stats_avg - avg[k]);

        // insert in-between for all lists
        itm.insert(itm.begin() + k + 1, setup.inv_T);
        avg.insert(avg.begin() + k + 1, stats_avg);
        var.insert(var.begin() + k + 1, stats_var);

        // if right interval is chosen, index points to inserted as lower bound
        k += (is_right);
    }
}

void MCMC::power_post_fractions_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
{
    for (setup.rung = 1; setup.rung <= setup.n_rungs; ++setup.rung)
    {
        // inverse tempearture = (rung / #rungs)^5
        setup.inv_T = std::pow(static_cast<double>(setup.rung) / static_cast<double>(setup.n_rungs), 5.);
        run_monte_carlo_();

        const auto [stats_avg, stats_var] = compute_stats(tot_log_lik.smpl_values);
        itm.push_back(setup.inv_T);
        avg.push_back(stats_avg);
        var.push_back(stats_var);
    }
}
*/
