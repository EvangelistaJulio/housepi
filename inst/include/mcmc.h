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
	};

	/// Engine running a built-in or custom model: evaluation, inference and simulation.
	class MCMC : public MCMCBase
	{
	public:
		std::unique_ptr<ModelBase> model; /// Transmission model.

		//////  Constructor  //////

		explicit MCMC(std::string_view model_name);

		/// Detach a shared data store before mutation (copy-on-write).
		void detach_data_store();

		//////  Evaluation  //////

		/// Initialize the model and evaluate the complete-data log-likelihood, joint over latents.
		/// - sum of observation, transmission and latent terms;
		/// - missing `hhid`: every household.
		double evaluate_log_lik(std::string_view hhid);

	private:
		/// Assert definitions, then initialize for evaluation.
		void init_for_evaluation_();

		/// Compute the total log-posterior, storing every old value.
		/// - result: `tot_obs_log_lik`, `tot_transm_log_lik`, `tot_latent_log_lik`;
		/// - `with_latent_ll`, `with_obs_ll`: re-evaluate the latent and observation terms.
		void commit_total_log_posterior_(bool with_cache,
										 bool with_latent_ll = true,
										 bool with_obs_ll = true);
		/// Commit the posterior contribution of `house`.
		LogLikAccumulator commit_house_log_posterior_(const Household &house,
													  bool with_cache,
													  bool with_latent_ll,
													  bool with_obs_ll);

		/// Compute the total log-posterior at initialization, guarded.
		/// - non-finite contribution: recovered in place, or error.
		void init_guarded_commit_total_();
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

	public:
		//////  Metropolis-Hastings  //////

		/// Run inference with the Metropolis-Hastings algorithm.
		/// @return Seconds spent sampling.
		double run_inference();

	private:
		/// Run every chain: initialization, adaptation, then posterior sampling.
		void run_monte_carlo_();
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
}
