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

#include "tools/nan.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  *** Total Log-Likelihood : struct ***  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

struct TotalLogLike
{
public:
    std::vector<double> smpl_values = {}; /// vector of total log-likelihood samples

    void recover() noexcept { value_ = buffer_; }
    double logratio() const noexcept { return value_ - buffer_; }

    /// Set new total log-likelihood value and buffer old value.
    void set(double new_ll) noexcept
    {
        buffer_ = value_;
        value_ = new_ll;
    }
    /// Increment total log-likelohood by a delta, buffer old value.
    void increment(double delta_ll) noexcept { set(value_ + delta_ll); }

    void init_draws(unsigned long sample_size)
    {
        count_ = 0;
        smpl_values.clear();
        smpl_values.assign(sample_size, hsp::na::DOUBLE);
    }
    void store_draws() { smpl_values[count_++] = value_; }
    unsigned long count() const noexcept { return count_; }

private:
    double value_ = 0.;       /// current total log-likelihood
    double buffer_ = 0.;      /// buffer total log-likelihood
    unsigned long count_ = 0; /// count of store samples
};

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

    TotalLogLike tot_log_lik;

    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    MCMCImpl(
        const std::vector<std::string> &data_wide_keys,
        const std::vector<std::string> &hh_size_t_keys,
        const std::vector<std::string> &hh_double_keys,
        const std::vector<std::string> &hh_matrix_keys,
        const std::vector<std::string> &data_long_keys,
        const std::vector<std::string> &param_keys,
        const std::vector<std::string> &group_rate_keys,
        const std::vector<std::string> &augm_keys,
        std::function<void()> func_init_inference,
        std::function<double(Household &)> func_house_log_lik,
        std::function<void(Household &, Individual &)> func_update_cache);
    ~MCMCImpl();

    ///////////////////////////////////
    //////  METROPOLIS-HASTINGS  //////
    ///////////////////////////////////

    void init_monte_carlo();
    void init_evaluate();
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

private:
    /////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////  MCMC : Impl : Callbacks Virtual Methods  ////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////

    /// UPDATE: this is not efficient.
    std::function<void()> callback_init_inference_;
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
