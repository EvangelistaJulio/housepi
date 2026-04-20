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
		/// UPDATE: can I transition to template now?
		std::unique_ptr<ModelBase> model;

		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		explicit MCMC(std::string_view model_name);

		/// Detach shared DataStore before mutation. Used for copy-on-write logic.
		void detach_data_store();

		/// Initialize frequency equivalence classes.
		/// NOTE: assume unique if augmented definitions required, but not imported yet
		void init_frequency_equivalence_classes();

		//////////////////////////
		//////  EVALUATION  //////
		//////////////////////////

		/// Initialize model and evaluate log-likelihood.
		double evaluate_log_lik(std::string_view house_id);

		///////////////////////////////////
		//////  METROPOLIS-HASTINGS  //////
		///////////////////////////////////

		/// Performs the inference using Metropolis-Hastings.
		void run_inference();

		/* --FUTURE-LOG-EVIDENCE--
		/// Trapezoidal rule to estimate logarithm of marginal likelihood (evidence) with power posteriors.
		///
		/// Ref
		/// -----
		/// Friel, N., Hurn, M. & Wyse, J.
		/// Improving power posterior estimation of statistical evidence.
		/// Stat Comput 24, 709–723 (2014).
		std::vector<std::vector<double>> power_posteriors(std::string_view scheme, unsigned long N);
		*/

	private:
		//////////////////////////
		//////  EVALUATION  //////
		//////////////////////////

		/// Assert setup definitions, then proceed with initialization for evaluation.
		void init_evaluate_(bool is_reset);

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

		/// Compute new log-likelihood, store old value. Same for households.
		double commit_total_log_lik_(bool is_cache_depen);
		
		///////////////////////////////////
		//////  METROPOLIS-HASTINGS  //////
		///////////////////////////////////

		/// Markov chain core, used for inference and power posteriors.
		/// NOTE: sweep is the number of previously completed thin sweeps.
		void run_monte_carlo_();

		/// Update and return contribution of given household to full logliklehood.
		double atomic_commit_house_log_lik_(Household &house, bool is_cache_depen);

		/// Move all model parameters and augmented variables.
		void move_all_();

		/// Attempt to move a single parameter for the given inverse temperature.
		bool move_single_parameter_(size_t p, bool is_cache_depen);

		/// Move augmented variables of all infected individuals.
		void move_all_augm_();

		/// Move augmented variables of infectious members inside household.
		/// The move is performed one member at a time, one augmented variable at a time.
		double atomic_move_augm_single_house_(size_t H);

		/// Move a single augmented variable.
		double metropolis_augm_single_move_(Household &house, Individual &member, size_t a);

		/* --FUTURE-LOG-EVIDENCE--
		/// Power posteriors at 0 inverse temperature.
		void power_post_zero_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
		/// Power posteriors with power fractions scheme for the thermodynamic ladder.
		void power_post_fractions_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
		void power_post_adaptive_(std::vector<double> &itm, std::vector<double> &avg, std::vector<double> &var);
		*/
	};
}
