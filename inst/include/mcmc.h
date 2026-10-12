/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "tools/logproba.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "mcmc/base.h"
#include "mcmc/setup.h"
#include "model.h"
#include "data/wide.h"
#include "data/long.h"
#include "parameter.h"
#include "risks.h"
#include "infec_profile.h"
#include "augmented.h"

//////////////////////////////  MCMC  //////////////////////////////

namespace housepi
{
	class Individual;
	class Household;

	/// Maximum resamples per member, at initialization.
	static constexpr size_t N_INIT_ATTEMPTS_ = 100;

	/// Sum log-likelihood terms over households.
	struct LogLikAccumulator
	{
		double obs_ll = 0.;	   /// Observation.
		double transm_ll = 0.; /// Transmission.
		double latent_ll = 0.; /// Latent model.

		LogLikAccumulator &operator+=(const LogLikAccumulator &other) noexcept
		{
			obs_ll += other.obs_ll;
			transm_ll += other.transm_ll;
			latent_ll += other.latent_ll;
			return *this;
		}
		/// Return the terms multiplied by `w`.
		LogLikAccumulator scaled(double w) const noexcept { return {w * obs_ll, w * transm_ll, w * latent_ll}; }
	};

	/// Engine running a built-in or custom model: evaluation, inference and simulation.
	class MCMC : public MCMCBase
	{
	public:
		std::unique_ptr<ModelBase> model; /// Transmission model.

		/// Create the engine of built-in or custom model `model_name`.
		explicit MCMC(std::string_view model_name);

		/// Detach a shared data store before mutation (copy-on-write).
		void detach_data_store();

		/// Initialize the model and evaluate the complete-data log-likelihood, joint over latents.
		/// - sum of observation, transmission and latent terms;
		/// - missing `hh_id`: every household, flagged evaluated.
		double evaluate_log_lik(std::string_view hh_id);

		/// Assert definitions and allocate the posterior, before sample().
		/// @note Call on the main thread: it may raise warnings.
		void prepare_inference();
		/// Run every chain: initialization, adaptation, then posterior sampling.
		/// - completion: every household flagged evaluated, at the last state.
		/// @return Seconds spent sampling.
		/// @note Call after prepare_inference(). Off the main thread, it calls no R.
		double sample();

	private:
		/// Create the engine of model `model_ptr`.
		explicit MCMC(std::unique_ptr<ModelBase> model_ptr);

		/// Assert definitions, then initialize for evaluation.
		void init_for_evaluation_();

		/// Compute the total log-posterior, storing every old value.
		/// - result: `tot_obs_log_lik`, `tot_transm_log_lik`, `tot_latent_log_lik`;
		/// - `with_latent_ll`, `with_obs_ll`: re-evaluate the latent and observation terms.
		void commit_total_log_posterior_(bool with_cache,
										 bool with_latent_ll = true,
										 bool with_obs_ll = true);
		/// Store the totals of `acc`.
		/// - transmission: always;
		/// - observation, latent model: when `do_obs`, `do_latent`.
		void set_total_log_lik_(const LogLikAccumulator &acc,
								bool do_obs,
								bool do_latent);
		/// Commit the posterior contribution of `house`.
		LogLikAccumulator commit_house_log_posterior_(const Household &house,
													  bool with_cache,
													  bool with_latent_ll,
													  bool with_obs_ll);

		/// Compute the guarded contribution of `house`.
		LogLikAccumulator init_guarded_commit_house_(Household &house);

		/// Compute and store the latent-model log-ratio of a move of `a`.
		double commit_latent_log_lik_ratio(const Individual &indiv, size_t a)
		{
			return augm.commit_latent_log_lik_ratio(indiv, a);
		}

		/// Compute and store the log-ratio of observations parameterized by `a`.
		double commit_obs_log_lik_ratio(const Individual &indiv, size_t a)
		{
			return augm.commit_obs_log_lik_ratio(indiv, a);
		}

		/// Compute and store the transmission log-likelihood of `house`.
		double commit_transm_log_lik_(const Household &house)
		{
			double output = model->compute_house_transm_log_lik(house);
			data.state.transm_log_lik(house).set(output);
			return output;
		}

		/// Compute and store the transmission log-likelihood of `house`.
		/// @return Log-ratio against the stored value.
		double commit_transm_log_lik_ratio_(const Household &house)
		{
			double output = -data.state.transm_log_lik(house).value();
			output += commit_transm_log_lik_(house);
			return output;
		}

		/// Initialize a chain: restart every variable, then commit the guarded total.
		/// - infeasible start: redraws, then error.
		void init_chain_()
		{
			init_monte_carlo_chain_();
			init_guarded_commit_total_();
		}
		/// Run chain `chain_id`: adaptation, then posterior sampling.
		/// @note Call after init_chain_(). It calls no R and raises no warning.
		void run_chain_(unsigned long chain_id)
		{
			run_mc_finite_adaptation_();
			run_mc_posterior_sampling_(chain_id);
		}

		/// Compute the total log-posterior at initialization, guarded.
		/// - non-finite contribution: recovered in place, or error.
		void init_guarded_commit_total_();
		/// Run the finite adaptation iterations.
		void run_mc_finite_adaptation_();
		/// Run posterior sampling iterations, burn-in included.
		/// - draws: stored in the column block of chain `chain_id`.
		void run_mc_posterior_sampling_(unsigned long chain_id);

		/// Move every parameter and latent.
		/// - `Adapt`: moves credited to the proposals, read by adaptation only.
		template <bool Adapt>
		void move_all_();
		/// Attempt a move of parameter `p`.
		template <bool Adapt>
		void move_single_parameter_(size_t p, bool with_cache);
		/// Move the latents of every infected member.
		template <bool Adapt>
		void move_all_latent_();
		/// Move the latents of infected members of household `H`, one member and one latent at a time.
		/// @return Accepted deltas of transmission and latent-model log-likelihoods.
		template <bool Adapt>
		LogLikAccumulator move_latent_single_house_(size_t H);
		/// Move latent `a` of member `indiv`.
		/// @return Accepted deltas of transmission and latent-model log-likelihoods.
		template <bool Adapt>
		LogLikAccumulator move_latent_single_variable_(const Household &house,
													   const Individual &indiv,
													   size_t a);
	};

	/// Outcome of run_inference().
	struct InferenceRun
	{
		std::vector<double> seconds; /// Seconds spent sampling, NA where the engine did not complete.
		parallel::Failure failure;	 /// Error that stopped the run, with the index of its engine.
	};

	/// Run inference on every engine, `workers` of them at once.
	/// Component-wise Metropolis-Hastings algorithm.
	/// Based on @cite metropolis_equation_1953 and @cite hastings_monte_1970
	///
	/// - preparation: in order, on the calling thread;
	/// - sampling: one engine per unit, in an arena of the largest `inner` among engines;
	/// - sampling on a new thread, `monitor` called on the calling thread every 250 ms;
	/// - error while preparing: nothing sampled;
	/// - error while sampling, or `monitor` throwing: every engine stops at its next iteration.
	/// @param monitor Callback taking the sum over engines of the fraction of iterations done.
	/// @note An exception thrown by `monitor` escapes once every engine has stopped.
	InferenceRun run_inference(const std::vector<MCMC *> &engine_list, size_t workers,
							   const std::function<void(double)> &monitor);
}
