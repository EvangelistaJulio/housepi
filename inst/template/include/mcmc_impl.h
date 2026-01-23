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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** MCMC : Impl : class ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::atomic
{
    struct Error;
}

class Individual;
class Household;
struct Setup;
struct DataWide;
struct DataLong;
struct Parameter;
struct GroupRate;
struct InfecProfile;
struct Augmented;

class MCMCImpl
{
public:
    std::unique_ptr<hsp::atomic::Error> atomic_error;
    std::shared_ptr<Setup> setup;
    std::shared_ptr<DataWide> data_wide;
    std::shared_ptr<DataLong> data_long;
    std::shared_ptr<Parameter> param;
    std::shared_ptr<GroupRate> group_rate;
    std::shared_ptr<InfecProfile> infec_profile;
    std::shared_ptr<Augmented> augm;

	///////////////////////////
	//////  CONSTRUCTOR  //////
	///////////////////////////

    MCMCImpl(std::function<void()> func_init_inference,
             std::function<void()> func_init_simulation,
             std::function<double(Household &)> func_house_log_lik,
             std::function<void(Household &)> func_simulate_house,
             std::function<void(Household &, Individual &)> func_update_cache);
    ~MCMCImpl();

	///////////////////////////////////
	//////  METROPOLIS-HASTINGS  //////
	///////////////////////////////////

    unsigned long smpl_count = 0;                /// count of store samples
    std::vector<double> smpl_total_log_lik = {}; /// vector of total log-likelihood samples

    void init_monte_carlo();
    void init_evaluate();
    void init_model_for_inference();
    void run_monte_carlo();
    std::vector<std::vector<double>> power_posteriors(std::string_view scheme, unsigned long N);
    double evaluate_single_house_log_lik(std::string_view house_id);
    double update_total_log_lik(bool is_cache_depen = true);

	//////////////////////////
	//////  SIMULATION  //////
	//////////////////////////

    /// Auxiliary method to initialize simulations.
    /// UPDATE: once I allow for multi-treading and a flag other than n_chains = 0, I can drop this.
    void helper_init_simulation();

    /// Simulate a single household outbreak in place.
    void simulate_single_house_in_place(Household &house);

private:
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  MCMC : Impl : Callbacks Virtual Methods  ////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    std::function<void()> callback_init_inference_;
    std::function<void()> callback_init_simulation_;
    std::function<double(Household &)> callback_house_log_lik_;
    std::function<void(Household &)> callback_simulate_house_;
    std::function<void(Household &, Individual &)> callback_update_cache_;

    void init_model_for_simulation_() { callback_init_simulation_(); }
    void update_cache_dependecies_(Household &h, Individual &i) { callback_update_cache_(h, i); }
    double materialize_house_obs_log_lik_(Household &h);

    /////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////  MCMC : Impl : Metropolis-Hastings [private]  //////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    double total_log_lik_ = 0.; /// current total log-likelihood
    double total_buffer_ = 0.;  /// buffer total log-likelihood

    /// Set new total log-likelihood value and buffer old value.
    void set_total_log_lik_(double new_ll)
    {
        total_buffer_ = total_log_lik_;
        total_log_lik_ = new_ll;
    }

    /// Increment total log-likelohood by a delta, buffer old value.
    void increment_total_log_lik_(double delta) { set_total_log_lik_(total_log_lik_ + delta); }

    void init_move_counts_();

    double atomic_materialize_house_log_lik_(Household &house, bool is_cache_depen = true);

    void move_all_();
    void adapt_proposal_of_sampled_params_(unsigned long step);

    void init_posterior_sample_();
    void store_posterior_sample_();

    //----- Metropolis - parameter moves -----//
    bool move_single_parameter_(size_t p, bool is_cache_depen);
    bool propose_valid_param_move_(size_t p);
    bool accept_valid_param_move_(size_t p);
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
