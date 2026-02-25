/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

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

#include "tools/types.h"

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** TotalLogLike : class ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

struct TotalLogLike : public hsp::LogLike
{
public:
	std::vector<double> smpl_values = {}; /// total log-likelihood samples

	void init_draws(unsigned long sample_size)
	{
		count_ = 0;
		smpl_values.clear();
		smpl_values.assign(sample_size, hsp::na::DOUBLE);
	}
	void store_draws() { smpl_values[count_++] = value_; }
	unsigned long count() const noexcept { return count_; }

private:
	unsigned long count_ = 0; /// count of stored samples
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  *** MCMC : class ***  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class AugmContainer;
class Individual;
class Household;

class MCMC
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

	TotalLogLike tot_log_lik;

	///////////////////////////
	//////  CONSTRUCTOR  //////
	///////////////////////////

	MCMC(const std::vector<std::string> &data_wide_keys,
		 const std::vector<std::string> &hh_size_keys,
		 const std::vector<std::string> &hh_size_vec_keys,
		 const std::vector<std::string> &hh_double_vec_keys,
		 const std::vector<std::string> &hh_double_mtx_keys,
		 const std::vector<std::string> &data_long_keys,
		 const std::vector<std::string> &param_keys,
		 const std::vector<std::string> &group_rate_keys,
		 const std::vector<std::string> &augm_keys);
	virtual ~MCMC();

	/////////////////////
	//////  MODEL  //////
	/////////////////////

	virtual void init_advanced_setup() {}

	/// Assign the contact patterns between a couple of individuals.
	/// UPDATE: this will be dropped in favor of table definition.
	virtual std::string init_contact_rates(Household &, const Individual &, const Individual &) { return ""; }

	/// Update cached augmeted dependecies for whole household.
	virtual void update_house_cache_dependecies(Household &) {}
	/// Update cached augmeted dependecies for single household member.
	virtual void update_member_cache_dependecies(Household &, Individual &) {}
	/// Compute the log-likelihood for a single household.
	virtual double compute_house_obs_log_lik(Household &) = 0;

	/// Inference specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_evaluation() {}
	/// Inference specific initialization for each household.
	/// NOTE: it's run right after 'init_model_for_evaluation'.
	virtual void init_house_for_evaluation(Household &) {}

	/// Simulation specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_simulation() {}
	/// Simulation for a single household.
	virtual void simulate_single_house(Household &) = 0;

	//////////////////////
	//////  SETTER  //////
	//////////////////////

	void init_infec_profile_template(hsp::dict<std::string> x) { infec_profile.init_template(x); }
	/// User call. Initialize all augmented variables of given individual.
	void init_missing_augm(Individual &member)
	{
		augm.init_missing(member, false);
		augm.user_reset(setup.prng(), member);
	}

	/// User call. Request discrete parameter definition.
	void request_discrete_param(size_t p) { param.request_discrete(p); }
	/// User call. Request discrete augmented definition.
	void request_discrete_augm(size_t a) { augm.request_discrete(a); }
	/// User call. Request discrete infectivity profile definition.
	void request_discrete_infec_profile() { infec_profile.request_discrete(); }

	/// Auxiliary method to initialize simulations.
	/// UPDATE: once I allow for multi-treading and a flag other than n_chains = 0, I can drop this.
	void helper_init_simulation()
	{
		data_wide.assert_non_empty(); // check if data is empty
		setup.n_chains = 0;			  // use 0 chains as a state variable
		setup.n_threads = 1;		  // impose no parallelization
		setup.set_n_threads();		  // change accordingly
	}

	//////////////////////
	//////  GETTER  //////
	//////////////////////

	/// Sample a uniform-distributed variable between 0 and 1.
	double uniform() const { return setup.uniform(); }
	/// Sample a exponential-distributed variable between 0 and Inf.
	double exponential() const { return setup.exponential(); }

	/// Get parameter value.
	double par(size_t p) const { return param.at(p).value(); }

	/// Return maximum number of observed infected household members.
	size_t get_max_hh_cases() const noexcept { return data_wide.get_max_hh_cases(); }
	/// Throw an exception if at least a member as missing status.
	void assert_no_missing_status() const { data_wide.assert_no_missing_status(); }

	double infec_profile_inst(const AugmContainer &i, double t) const { return infec_profile.inst(&i, t); }
	double infec_profile_cumul(const AugmContainer &i, double t) const { return infec_profile.cumul(&i, t); }

	/// Return cached size_t value.
	size_t &cache_size(Household &house, size_t e) { return data_wide.cache_size(house, e); }
	size_t cache_size(const Household &house, size_t e) const { return data_wide.cache_size(house, e); }

	/// Return cached size_t vector.
	std::vector<size_t> &cache_size_vector(Household &house, size_t e) { return data_wide.cache_size_vector(house, e); }
	const std::vector<size_t> &cache_size_vector(const Household &house, size_t e) const { return data_wide.cache_size_vector(house, e); }

	/// Return cached double vector.
	std::vector<double> &cache_double_vector(Household &house, size_t e) { return data_wide.cache_double_vector(house, e); }
	const std::vector<double> &cache_double_vector(const Household &house, size_t e) const { return data_wide.cache_double_vector(house, e); }

	/// Return cached double matrix.
	std::vector<std::vector<double>> &cache_double_matrix(Household &house, size_t e) { return data_wide.cache_double_matrix(house, e); }
	const std::vector<std::vector<double>> &cache_double_matrix(const Household &house, size_t e) const { return data_wide.cache_double_matrix(house, e); }

	//////////////////////////
	//////  EVALUATION  //////
	//////////////////////////

	/// Initialize model and evaluate log-likelihood.
	double evaluate_log_lik(std::string_view house_id = "");

	///////////////////////////////////
	//////  METROPOLIS-HASTINGS  //////
	///////////////////////////////////

	/// Performs the inference using Metropolis-Hastings.
	void run_inference();
	
	/// Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
	///
	/// Ref
	/// -----
	/// Friel, N., Hurn, M. & Wyse, J.
	/// Improving power posterior estimation of statistical evidence.
	/// Stat Comput 24, 709–723 (2014).
	std::vector<std::vector<double>> power_posteriors(std::string_view scheme, unsigned long N);

private:
	//////////////////////////
	//////  EVALUATION  //////
	//////////////////////////

	/// Assert setup definitions, then proceed with initialization for evaluation.
	void init_evaluate_(bool is_reset = true);
	/// Compute new log-likelihood, store old value in buffer. Same for households.
	double update_total_log_lik_(bool is_cache_depen = true);

	///////////////////////////////////
	//////  METROPOLIS-HASTINGS  //////
	///////////////////////////////////

	/// Assert setup definitions are ok to proceed with mcmc sampling.
	void init_monte_carlo_();
	/// Intialize a single Markov chain for Monte Carlo.
	void init_monte_carlo_chain_();
	/// Markov chain core, used for inference and power posteriors.
	/// NOTE: sweep is the number of previously completed thin sweeps.
	void run_monte_carlo_();

	/// Update and return contribution of given household to full logliklehood.
	double atomic_materialize_house_log_lik_(Household &house, bool is_cache_depen = true);

	/// Move all model parameters and augmented variables.
	void move_all_();
	/// Attempt to move a single parameter for the given inverse temperature.
	bool move_single_parameter_(size_t p, bool is_cache_depen);
	/// Restore the parameter value in the buffer.
	void recover_parameter_(size_t p, bool is_cache_depen);
	/// Move augmented variables of all infected individuals.
	void move_all_augm_();
	/// Move augmented variables of infectious members inside household.
	/// The move is performed one member at a time, one augmented variable at a time.
	double atomic_move_augm_single_house_(size_t H);
	/// Move a single augmented variable.
	double metropolis_augm_single_move_(Household &house, Individual &member, size_t a);

	/// Init posterior sampling.
	//// UPDATE: I would like to store also power posterior sample.
	/// NOTE: it will depend on chosen rungs in that case.
	void init_posterior_sample_();
	/// Store current Bayesian posterior sample.
	void store_posterior_sample_();

	/// Power posteriors at 0 inverse temperature.
	void power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
	/// Power posteriors with power fractions scheme for the thermodynamic ladder.
	void power_post_fractions_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
	void power_post_adaptive_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
};
