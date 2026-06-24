/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "tools/logproba.h"
#include "tools/nan.h"

#include "mcmc_base.h"
#include "model.h"
#include "setup.h"
#include "data_wide.h"
#include "data_long.h"
#include "parameter.h"
#include "group_rate.h"
#include "infec_profile.h"
#include "augmented.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  *** MCMC : class ***  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
	class Individual;
	class Household;

	/// Output of a power-posterior ladder run.
	struct PowerLadder
	{
		std::vector<double> inv_temp;			/// inverse-temperature ladder
		std::vector<double> mean_log_lik;		/// per-rung mean log-likelihood
		std::vector<double> var_log_lik;		/// per-rung variance log-likelihood
		std::vector<std::vector<double>> draws; /// total log-likelihood draw

		PowerLadder(size_t n_rungs); /// reserve all vectors

		void push(double inv_t, std::vector<double> row); /// compute stats + append
		void insert(size_t i, double inv_t, std::vector<double> row);
	};

	class MCMC : public MCMCBase
	{
	public:
		std::unique_ptr<ModelBase> model;

		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		explicit MCMC(std::string_view model_name);

		/// Detach shared DataStore before mutation. Used for copy-on-write logic.
		void detach_data_store();

		//////////////////////////
		//////  EVALUATION  //////
		//////////////////////////

		/// Initialize model and evaluate the scalar complete-data household
		/// log-likelihood (observation + transmission), joint over augmented variables.
		double evaluate_log_lik(std::string_view house_id);

	private:
		/// Assert setup definitions, then proceed with initialization for evaluation.
		void init_for_evaluation_();

		/// Compute total log-posterior. Store all old values, included for households.
		/// Result is left in tot_transm_log_lik / tot_augm_log_prior / tot_obs_log_lik.
		/// NOTE: with_augm_prior / with_obs skip the respective re-evaluation.
		void commit_total_log_posterior_(bool is_cache_depen,
										 bool with_augm_prior = true,
										 bool with_obs = true);
		/// Commit single household's posterior contribution.
		LogProbaAccumulator commit_house_log_posterior_(Household &house,
														bool is_cache_depen,
														bool with_augm_prior,
														bool with_obs);

		/// Guarded total log-posterior computation for initialization.
		/// Assert for finite initial value and mechanisms to recover in place.
		void init_guarded_commit_total_();
		/// Single household's guarded contribution.
		LogProbaAccumulator init_guarded_commit_house_(Household &house);

		/// Compute and store log-ratio of prior of single augmented variable.
		double commit_augm_log_prior_ratio(Individual &member, size_t a)
		{
			return augm.commit_augm_log_prior_ratio(member, a);
		}

		/// Compute and store log-ratio of observation nodes parameterized by 'a'.
		double commit_augm_log_obs_ratio(Individual &member, size_t a)
		{
			return augm.commit_obs_log_lik_ratio(member, a);
		}

		/// Compute and set value of observation log-likelihood. Return log-like.
		double commit_transm_log_lik_(Household &house)
		{
			double output = model->compute_house_transm_log_lik(house);
			data_wide().transm_log_lik(house).set(output);
			return output;
		}

		/// Compute and set value of observation log-likelihood. Return log-like ratio.
		double commit_transm_log_lik_ratio_(Household &house)
		{
			double output = -data_wide().transm_log_lik(house).value();
			output += commit_transm_log_lik_(house);
			return output;
		}

	public:
		///////////////////////////////////
		//////  METROPOLIS-HASTINGS  //////
		///////////////////////////////////

		/// Performs the inference using Metropolis-Hastings algorithm.
		void run_inference();

		/// Run power-posterior thermodynamic ladder.
		/// Based on @cite friel_improving_2014
		PowerLadder run_power_posteriors(std::string_view scheme, unsigned long N);

	private:
		/// Markov chain core, used for inference and power posteriors.
		/// NOTE: sweep is the number of previously completed thin sweeps.
		void run_monte_carlo_();
		/// Finite adaptative iterations.
		void run_mc_finite_adaptation_();
		// Run MC chain posterior sampling iterations, with optional burn-in.
		void run_mc_posterior_sampling_();

		/// Move all model parameters and augmented variables.
		void move_all_();
		/// Attempt to move a single parameter.
		void move_single_parameter_(size_t p, bool is_cache_depen);
		/// Move augmented variables of all infected individuals.
		void move_all_augm_();
		/// Move augmented variables of infectious members inside household.
		/// The move is performed one member at a time, one augmented variable at a time.
		/// Returns the accepted transmission log-likelihood and augmented log-prior deltas.
		LogProbaAccumulator move_augm_single_house_(size_t H);
		/// Move a single augmented variable.
		/// Returns the accepted transmission log-likelihood and augmented log-prior deltas.
		LogProbaAccumulator move_augm_single_variable_(Household &house, Individual &member, size_t a);

		/// Power posteriors at 0 inverse temperature.
		void power_post_zero_(PowerLadder &ladder);
		/// Power posteriors with power fractions scheme for the thermodynamic ladder.
		void power_post_fractions_(PowerLadder &ladder);
		void power_post_adaptive_(PowerLadder &ladder);
	};
}
