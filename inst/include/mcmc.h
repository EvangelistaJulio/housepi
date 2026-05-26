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

#include "tools/loglike.h"
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

	class MCMC : public MCMCBase
	{
	public:
		/// UPDATE: can I drop pointer and just transition to template?
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

		/// Initialize model and evaluate log-likelihood.
		double evaluate_log_lik(std::string_view house_id);

	private:
		/// Assert setup definitions, then proceed with initialization for evaluation.
		void init_evaluate_(bool is_reset);

		/// Compute new log-likelihood, store old value. Same for households.
		double commit_total_log_lik_(bool is_cache_depen);
		/// Update and return contribution of given household to full logliklehood.
		double commit_house_log_lik_(Household &house, bool is_cache_depen);
		/// Compute and set value of observation log-likelihood for given household.
		double commit_house_obs_log_lik_(Household &house)
		{
			double output = model->compute_house_obs_log_lik(house);
			house.obs_log_lik_.set(output);
			return output;
		}
		/// Compute and set value of observation log-likelihood for given household.
		double commit_house_obs_log_ratio_(Household &house)
		{
			double output = -house.obs_log_lik_.value();
			output += commit_house_obs_log_lik_(house);
			return output;
		}

	public:
		///////////////////////////////////
		//////  METROPOLIS-HASTINGS  //////
		///////////////////////////////////

		/// Performs the inference using Metropolis-Hastings.
		void run_inference();

		/* --FUTURE-LOG-EVIDENCE--
		/// Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
		/// Based on @cite friel_improving_2014
		std::vector<std::vector<double>> power_posteriors(std::string_view scheme, unsigned long N);
		*/

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
		/// Attempt to move a single parameter for the given inverse temperature.
		void move_single_parameter_(size_t p, bool is_cache_depen);
		/// Move augmented variables of all infected individuals.
		void move_all_augm_();
		/// Move augmented variables of infectious members inside household.
		/// The move is performed one member at a time, one augmented variable at a time.
		double move_augm_single_house_(size_t H);
		/// Move a single augmented variable.
		double move_augm_single_variable_(Household &house, Individual &member, size_t a);

		/* --FUTURE-LOG-EVIDENCE--
		/// Power posteriors at 0 inverse temperature.
		void power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
		/// Power posteriors with power fractions scheme for the thermodynamic ladder.
		void power_post_fractions_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
		void power_post_adaptive_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
		*/
	};
}
