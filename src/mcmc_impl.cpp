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
#include "registry.h"
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

namespace hsp
{
    MCMCBase::MCMCBase() : pDataStore(std::make_shared<DataStore>()), augm(this) {}
    MCMCBase::~MCMCBase() = default;

    void MCMCBase::set_data_store(std::shared_ptr<DataStore> ds)
    {
        pDataStore = std::move(ds);
        DataWide *dw = &pDataStore->data_wide;
        param.bind_to_data_wide(dw->pKeys.get());
    }

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
        data_wide().init_registry(data_wide_keys,
                                  hh_size_keys,
                                  hh_size_vec_keys,
                                  hh_double_vec_keys,
                                  hh_double_mtx_keys);
        data_long().init_registry(data_long_keys);
        param.init_registry(param_keys);
        group_rate.init_registry(group_rate_keys);
        augm.init_registry(augm_keys);

        /// UPDATE: waiting for refactoring of read_attribute
        param.bind_to_data_wide(data_wide().pKeys.get());
        param.bind_to_group_rate(group_rate.pKeys.get());
        param.bind_to_augm(augm.pKeys.get());
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
        // recover old value of parameter
        param.at(p).recover();
        // recover old value of full log-like
        tot_log_lik.recover();

        // recover old log-likelihood values for households
        for (size_t R = 0; R < data_wide().restricted_size(); ++R)
        {
            Household &house = data_wide().restricted(R);
            house.obs_log_lik_.recover();

            if (augm.is_any_sample())
                for (size_t i : house.infected_hhindexes())
                    house.member(i).paugm_->log_lik.recover();
        }

        if (is_cache_depen)
        {
            data_wide().recover_old_cached();
            infec_profile.recover();
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  *** MCMC : class ***  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    MCMC::MCMC(std::string_view model_name)
        : model(ModelRegistry::instance().create(model_name))
    {
        init_registry(model->get_DataKey(),
                      model->get_HHsizeKey(),
                      model->get_HHsizevecKey(),
                      model->get_HHdoublevecKey(),
                      model->get_HHdoublemtxKey(),
                      model->get_LongKey(),
                      model->get_ParamKey(),
                      model->get_GroupRateKey(),
                      model->get_AugmKey());
        model->bind_solver(this);
        model->init_advanced_setup();
    }

    void MCMC::detach_data_store()
    {
        if (pDataStore.use_count() <= 1)
            return;
        auto new_ds = std::make_shared<DataStore>();
        new_ds->data_wide.init_registry(
            model->get_DataKey(),
            model->get_HHsizeKey(),
            model->get_HHsizevecKey(),
            model->get_HHdoublevecKey(),
            model->get_HHdoublemtxKey());
        new_ds->data_long.init_registry(model->get_LongKey());
        set_data_store(std::move(new_ds));
    }

    //////////////////////////
    //////  EVALUATION  //////
    //////////////////////////

    double MCMC::evaluate_log_lik(std::string_view house_id)
    {
        setup.n_chains = 0;   // use 0 chains as a state variable
        init_evaluate_(true); // initialize for evaluation

        // mode: single house evaluation
        if (!na::is(house_id))
        {
            infec_profile.update(); // eventually recompute infectivity profile
            double output = atomic_commit_house_log_lik_(data_wide().at(house_id), true);
            atomic_error.throw_exception(); // catch and throw if previous failed
            return output;
        }

        // mode: total evaluation
        return commit_total_log_lik_(true);
    }

    ///////////////////////////////////
    //////  METROPOLIS-HASTINGS  //////
    ///////////////////////////////////

    void MCMC::run_inference()
    {
        init_monte_carlo_();               // initialize mcmc sampling
        init_evaluate_(false);             // initialize for evaluation
        setup.prog_open_for_mcmc_sample(); // open progress bar
        run_monte_carlo_();                // run Monte Carlo
        setup.prog_close();                // close progress bar
    }

    //////////////////////////
    //////  EVALUATION  //////
    //////////////////////////

    void MCMC::init_evaluate_(bool is_reset)
    {
        data_wide().assert_non_empty();                      // check if data is empty
        data_wide().assert_not_all_missing_status();         // warn if all missing status
        augm.init();                                         // initialize individual augmented relations
        data_wide().init_multiplicity(augm.is_any_sample()); // multiplicity

        if (is_reset)
        {
            param.user_reset(setup.prng()); // initialize non constant parameters
            augm.user_reset(setup.prng());  // initialize non constant augmented
        }

        model->init_model_for_evaluation(); // user initialization of model for inference
        data_wide().allocate_all_cache();     // allocate all cache dependencies

        for (size_t R = 0; R < data_wide().restricted_size(); ++R)
            model->init_house_for_evaluation(data_wide().restricted(R));
    }

    double MCMC::commit_total_log_lik_(bool is_cache_depen)
    {
        double new_ll = 0.0;

        // save infectivity profile if lookup tables are active
        if (is_cache_depen)
        {
            data_wide().save_all_cached();
            infec_profile.update();
        }

#pragma omp parallel for reduction(+ : new_ll)
        for (size_t R = 0; R < data_wide().restricted_size(); ++R)
        {
            if (atomic_error.is_pending())
                continue;

            Household &house = data_wide().restricted(R);
            new_ll += data_wide().multiplicity(R) *
                      atomic_commit_house_log_lik_(house, is_cache_depen);
        }
        atomic_error.throw_exception();
        tot_log_lik.set(new_ll);

        return new_ll;
    }

    ///////////////////////////////////
    //////  METROPOLIS-HASTINGS  //////
    ///////////////////////////////////

    void MCMC::run_monte_carlo_()
    {
        posterior.init(*this); // initialize posterior sample

        // run multiple MCMC chains&
        for (setup.chain = 1; setup.chain <= setup.n_chains; ++setup.chain)
        {
            try
            {
                init_monte_carlo_chain_();
                // NOTE: init all cached with flag `true`
                commit_total_log_lik_(true);

                // assert that initial log-likelihood value is finite
                if (tot_log_lik.is_inf())
                    data_wide().assert_finite_log_lik();
            }
            catch (const std::exception &e)
            {
                throw std::runtime_error(
                    "While initializing Markov chain:\n" + std::string(e.what()));
            }

            // Finite adaptative iterations
            for (unsigned long i = 0; i < setup.n_adapt; i += setup.n_thin, setup.prog_update())
                for (unsigned long t = 1; t <= setup.n_thin; ++t)
                {
                    move_all_();
                    param.adapt_proposal_of_sampled(i + t, setup.n_adapt);
                }

            // Burn-in iterations
            for (unsigned long i = 0; i < setup.n_burn; i += setup.n_thin, setup.prog_update())
                for (unsigned long t = 1; t <= setup.n_thin; ++t)
                    move_all_();

            // Posterior sampling
            for (unsigned long i = setup.n_burn; i < setup.n_steps; i += setup.n_thin, setup.prog_update())
            {
                for (unsigned long t = 1; t <= setup.n_thin; ++t)
                    move_all_();
                posterior.store(*this);
            }
        }
    }

    double MCMC::atomic_commit_house_log_lik_(Household &house, bool is_cache_depen)
    {
        double output = 0;
        try
        {
            // save all household cache dependencies
            if (is_cache_depen)
                model->update_house_cache_dependecies(house);

            // compute and set household observation log-likelihood
            output = commit_house_obs_log_lik_(house);

            // compute and set individual augmented log-likelihood
            if (augm.is_any_model())
                for (size_t i : house.infected_hhindexes())
                    output += augm.commit_augm_model_log_lik(house.member(i));
        }
        catch (const std::exception &e)
        {
            atomic_error.set("While computing log-likelihood for " +
                             house.err_name() + ":\n" + std::string(e.what()));
        }

        return output;
    }

    void MCMC::move_all_()
    {
        // move augmented variables to sample, if any
        if (augm.is_any_sample())
            move_all_augm_();

        // move all standard model parameters
        for (size_t p : param.list_adaptable_sample_standard)
            param.at(p).update_accept_rate(move_single_parameter_(p, false));
        for (size_t p : param.list_non_adapt_sample_standard)
            move_single_parameter_(p, false);

        // move parameters with cached dependecies
        for (size_t p : param.list_adaptable_sample_cached)
            param.at(p).update_accept_rate(move_single_parameter_(p, true));
        for (size_t p : param.list_non_adapt_sample_cached)
            move_single_parameter_(p, true);
    }

    bool MCMC::move_single_parameter_(size_t p, bool is_cache_depen)
    {
        // proceed only if the new value is consistent
        // otherwise automatically refuse move
        if (param.propose_valid_move(p, setup.prng()))
        {
            // if parameter is cache dependent, update cache
            commit_total_log_lik_(is_cache_depen);

            // perform Metropolis-Hastings to accept or not move
            const double r = setup.inv_T * tot_log_lik.log_ratio() +
                             param.at(p).get_move_log_ratio();

            if (setup.exponential() > -r)
                return true; // accepted move
            recover_parameter_(p, is_cache_depen);
        }
        return false; // refused move
    }

    void MCMC::move_all_augm_()
    {
        double delta_ll = 0.;

#pragma omp parallel for reduction(+ : delta_ll)
        for (size_t H = 0; H < data_wide().report_house.total; ++H)
        {
            if (atomic_error.is_pending())
                continue;

            delta_ll += atomic_move_augm_single_house_(H);
        }
        atomic_error.throw_exception();
        tot_log_lik.increment(delta_ll);
    }

    double MCMC::atomic_move_augm_single_house_(size_t H)
    {
        Household &house = data_wide().at(H);
        double delta_ll = 0.;

        try
        {
            // loop over augmented variables of infectious
            for (size_t i : house.infected_hhindexes())
            {
                Individual &member_i = house.member(i);
                for (size_t a : augm.list_sample(member_i))
                    delta_ll += metropolis_augm_single_move_(house, member_i, a);
            }
        }
        catch (const std::exception &e)
        {
            atomic_error.set("While sampling augmented variables for " +
                             house.err_name() + ":\n" + std::string(e.what()));
        }

        return delta_ll;
    }

    double MCMC::metropolis_augm_single_move_(Household &house, Individual &member, size_t a)
    {
        Variable &xaugm = member.get_augm_variable(a);

        xaugm.propose(setup.prng()); // propose new value and save old

        // proceed only if the new value is consistent
        if (xaugm.is_valid_move())
        {
            size_t i = member.hhindex();
            // update dependencies to perform computations
            data_wide().save_member_cache(house, i);
            model->update_member_cache_dependecies(house, member);

            // Metropolis-Hastings
            const double delta_ll = commit_house_obs_log_ratio_(house) +
                                    augm.commit_augm_model_log_ratio(member, a);

            if (setup.exponential() > xaugm.get_propose_log_ratio() - setup.inv_T * delta_ll)
                return delta_ll;

            // else recover cached values
            member.paugm_->log_lik.recover();
            house.obs_log_lik_.recover();
            data_wide().recover_member_old_cached(house, i);
        }
        // else recover old value
        xaugm.recover();

        return 0.;
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
                str::find_closest(scheme, options));

        // check number of rungs for the thermodynamic ladder
        if (n_rungs < 2)
            throw std::invalid_argument(
                "N. of rungs was not greater than 1.");

        setup.n_rungs = n_rungs;          // NOTE: set before progress bar
        double old_inv_T = setup.inv_T;   // save user setup
        init_monte_carlo_();              // initialize mcmc sampling
        init_evaluate_(false);            // initialize for evaluation
        setup.prog_open_for_power_post(); // print progress bar

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
        setup.inv_T = old_inv_T; // restore user setup

        return {itm, avg, var};
    }

    void MCMC::power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var)
    {
        // UPDATE: should I use previous inference for itm = 1 if available?
        // UPDATE: always under-estimate time if starting from itm = 0

        setup.rung = 0;
        setup.inv_T = 0.0;
        run_monte_carlo_();
        const auto [stats_avg, stats_var] = compute_stats(posterior.total_log_lik);

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
            const auto [stats_avg, stats_var] = compute_stats(posterior.total_log_lik);

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
            const auto [stats_avg, stats_var] = compute_stats(posterior.total_log_lik);

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

            const auto [stats_avg, stats_var] = compute_stats(posterior.total_log_lik);
            itm.push_back(setup.inv_T);
            avg.push_back(stats_avg);
            var.push_back(stats_var);
        }
    }
    */
}
