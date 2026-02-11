/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "tools/atomic.h"
#include "tools/loglike.h"
#include "tools/nan.h"

#include "setup.h"
#include "data_wide.h"
#include "data_long.h"
#include "parameter.h"
#include "group_rate.h"
#include "contact_rate.h"
#include "infec_profile.h"
#include "augmented.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** MCMC : Impl : class ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class AugmContainer;
class Individual;
class Household;

class MCMCImpl
{
public:
    hsp::atomic::Error atomic_error;

    Setup setup;
    DataWide data_wide;
    DataLong data_long;
    Parameter param;
    GroupRate group_rate;
    ContactRate contact_rate;
    InfecProfile infec_profile;
    Augmented augm;

    struct TotalLogLike : public hsp::LogLike
    {
    public:
        std::vector<double> smpl_values = {}; /// vector of total log-likelihood samples

        void init_draws(unsigned long sample_size)
        {
            count_ = 0;
            smpl_values.clear();
            smpl_values.assign(sample_size, hsp::na::DOUBLE);
        }
        void store_draws() { smpl_values[count_++] = value_; }
        unsigned long count() const noexcept { return count_; }

    private:
        unsigned long count_ = 0; /// count of store samples
    } tot_log_lik;

    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    MCMCImpl(
        const std::vector<std::string> &data_wide_keys,
        const std::vector<std::string> &hh_size_keys,
        const std::vector<std::string> &hh_size_vec_keys,
        const std::vector<std::string> &hh_double_vec_keys,
        const std::vector<std::string> &hh_double_mtx_keys,
        const std::vector<std::string> &data_long_keys,
        const std::vector<std::string> &param_keys,
        const std::vector<std::string> &group_rate_keys,
        const std::vector<std::string> &augm_keys,
        std::function<void()> func_init_inference,
        std::function<void(Household &)> func_init_house_for_infer,
        std::function<double(Household &)> func_house_log_lik,
        std::function<void(Household &, Individual &)> func_update_cache,
        std::function<std::string(Household &, const Individual &, const Individual &)> func_init_contact_rates);
    ~MCMCImpl();

    //////////////////////
    //////  GETTER  //////
    //////////////////////

    /// Sample a uniform-distributed variable between 0 and 1.
    double uniform() { return setup.uniform(); }
    /// Sample a exponential-distributed variable between 0 and Inf.
    double exponential() { return setup.exponential(); }

    double par(size_t p) const { return param.at(p).value(); }

    size_t get_max_hh_cases() const noexcept { return data_wide.get_max_hh_cases(); }
    void assert_no_missing_status() const { data_wide.assert_no_missing_status(); }

    double infec_profile_inst(const AugmContainer &i, double t) const { return infec_profile.inst(&i, t); }
    double infec_profile_cumul(const AugmContainer &i, double t) const { return infec_profile.cumul(&i, t); }

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// UPDATE: also contact rates should just be imported
    void build_contact_rates() { contact_rate.build(); }

    
    void init_infec_profile_template(hsp::dict<std::string> x) { infec_profile.init_template(x); }

    /// Initialize all augmented variables of given individual.
    void init_missing_augm(Individual &member);

    /// Request discrete parameter definition.
    void request_discrete_param(size_t p) { param.request_discrete(p); }
    /// Request discrete augmented definition.
    void request_discrete_augm(size_t a) { augm.request_discrete(a); }
    /// Request discrete infectivity profile definition.
    void request_discrete_infec_profile() { infec_profile.request_discrete(); }

    //////////////////////////
    //////  EVALUATION  //////
    //////////////////////////

    void user_init_for_evaluate();
    void init_evaluate();
    double update_total_log_lik(bool is_cache_depen = true);
    double evaluate_single_house_log_lik(std::string_view house_id);
    double evaluate_log_lik(std::string_view house_id = "");

    ///////////////////////////////////
    //////  METROPOLIS-HASTINGS  //////
    ///////////////////////////////////

    void init_monte_carlo();
    void run_monte_carlo();
    /// Performs the inference using Metropolis-Hastings.
    void run_inference();
    std::vector<std::vector<double>> power_posteriors(std::string_view scheme, unsigned long N);

    //////////////////////////
    //////  SIMULATION  //////
    //////////////////////////

    /// Auxiliary method to initialize simulations.
    /// UPDATE: once I allow for multi-treading and a flag other than n_chains = 0, I can drop this.
    void helper_init_simulation();

private:
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  MCMC : Impl : Callbacks Virtual Methods  ////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// UPDATE: this is not efficient.
    std::function<void()> callback_init_inference_;
    std::function<void(Household &)> callback_init_house_for_inference_;
    std::function<double(Household &)> callback_house_log_lik_;
    std::function<void(Household &, Individual &)> callback_update_cache_;

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  MCMC : Impl : Metropolis-Hastings [private]  //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// Intialize a single Markov chain for Monte Carlo.
    void init_monte_carlo_chain_();

    /// Update and return contribution of given household to full logliklehood.
    double atomic_materialize_house_log_lik_(Household &house, bool is_cache_depen = true);

    /// Move all model parameters and augmented variables.
    void move_all_();

    void init_posterior_sample_();
    void store_posterior_sample_();

    //----- Metropolis - parameter moves -----//
    bool move_single_parameter_(size_t p, bool is_cache_depen);
    void recover_parameter_(size_t p, bool is_cache_depen);

    //----- Metropolis : augmented moves -----//
    void move_all_augm_();
    double atomic_move_augm_single_house_(size_t H);
    double metropolis_augm_single_move_(Household &house, Individual &member, size_t a);

    //----- Metropolis : power posteriors -----//
    void power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
    void power_post_fractions_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
    void power_post_adaptive_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
};
