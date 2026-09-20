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
#include "group_risk.h"
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

	/// Max resamples per member.
	static constexpr size_t N_INIT_ATTEMPTS_ = 100;

	/// Output of a power-posterior ladder run.
	struct PowerLadder
	{
		std::vector<double> inv_temp;			/// inverse-temperature ladder
		std::vector<double> mean_log_lik;		/// per-rung mean log-likelihood
		std::vector<double> var_log_lik;		/// per-rung variance log-likelihood
		std::vector<std::vector<double>> draws; /// total log-likelihood draw
		std::string scheme;						/// scheme actually used to sample

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
		/// log-likelihood (observation + transmission), joint over latent variables.
		double evaluate_log_lik(std::string_view hhid);

	private:
		/// Assert setup definitions, then proceed with initialization for evaluation.
		void init_for_evaluation_();

		/// Compute total log-posterior. Store all old values.
		/// Result is left in tot_obs_log_lik, tot_transm_log_lik, tot_latent_log_lik.
		/// NOTE: with_latent_ll / with_obs_ll run respective re-evaluations.
		void commit_total_log_posterior_(bool with_cache,
										 bool with_latent_ll = true,
										 bool with_obs_ll = true);
		/// Commit single household's posterior contribution.
		LogLikAccumulator commit_house_log_posterior_(const Household &house,
													  bool with_cache,
													  bool with_latent_ll,
													  bool with_obs_ll);
		/// Recompute and store the total latent log-likelihood across all households.
		double commit_total_latent_log_lik_ratio_();
		/// Recover the total latent log-likelihood after a refused prior-only move.
		void recover_total_latent_log_lik_();

		/// Guarded total log-posterior computation for initialization.
		/// Assert for finite initial value and mechanisms to recover in place.
		void init_guarded_commit_total_();
		/// Single household's guarded contribution.
		LogLikAccumulator init_guarded_commit_house_(Household &house);

		/// Compute and store log-ratio of prior of single latent variable.
		double commit_latent_log_lik_ratio(const Individual &indiv, size_t a)
		{
			return augm.commit_latent_log_lik_ratio(indiv, a);
		}

		/// Compute and store log-ratio of observation nodes parameterized by 'a'.
		double commit_obs_log_lik_ratio(const Individual &indiv, size_t a)
		{
			return augm.commit_obs_log_lik_ratio(indiv, a);
		}

		/// Compute and set value of observation log-likelihood. Return log-like.
		double commit_transm_log_lik_(const Household &house)
		{
			double output = model->compute_house_transm_log_lik(house);
			data_wide().transm_log_lik(house).set(output);
			return output;
		}

		/// Compute and set value of observation log-likelihood. Return log-like ratio.
		double commit_transm_log_lik_ratio_(const Household &house)
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
		/// Beta0: chain at inverse temperature 0 (likelihood untempered, prior-only target).
		template <bool Beta0>
		void run_mc_finite_adaptation_();
		// Run MC chain posterior sampling iterations, with optional burn-in.
		// Draws are stored in the column block of chain `chain_id`.
		template <bool Beta0>
		void run_mc_posterior_sampling_(unsigned long chain_id);

		/// Move all model parameters and latent variables.
		/// Adapt: moves are credited to the proposals, which only adaptation reads.
		template <bool Adapt, bool Beta0>
		void move_all_();
		/// Attempt to move a single parameter.
		template <bool Adapt, bool Beta0>
		void move_single_parameter_(size_t p, bool with_cache);
		/// Move latent variables of all infected individuals.
		template <bool Adapt, bool Beta0>
		void move_all_latent_();
		/// Move latent variables of infectious members inside household.
		/// The move is performed one member at a time, one latent variable at a time.
		/// Returns the accepted transmission log-likelihood and latent log-model deltas.
		template <bool Adapt, bool Beta0>
		LogLikAccumulator move_latent_single_house_(size_t H);
		/// Move a single latent variable.
		/// Returns the accepted transmission log-likelihood and latent log-model deltas.
		template <bool Adapt, bool Beta0>
		LogLikAccumulator move_latent_single_variable_(const Household &house,
													 const Individual &indiv,
													 size_t a);

		/// Power posteriors at 0 inverse temperature.
		void power_post_zero_(PowerLadder &ladder);
		/// Power posteriors with power fractions scheme for the thermodynamic ladder.
		void power_post_fractions_(PowerLadder &ladder);
		void power_post_adaptive_(PowerLadder &ladder);
	};
}
