/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>
#include <string>

#include "tools/types.h"

#include "mcmc_base.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** ModelBase : class ***  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
	class Individual;
	class Household;

	class ModelBase
	{
	public:
		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		ModelBase() = default;
		virtual ~ModelBase() = default;

		void bind_solver(MCMCBase *mcmc) noexcept { mcmc_ = mcmc; }

		virtual std::vector<std::string> compile_flags() const { return {""}; }
		virtual const std::vector<std::string> get_DataKey() const = 0;
		virtual const std::vector<std::string> get_HHsizeKey() const = 0;
		virtual const std::vector<std::string> get_HHsizevecKey() const = 0;
		virtual const std::vector<std::string> get_HHdoublevecKey() const = 0;
		virtual const std::vector<std::string> get_HHdoublemtxKey() const = 0;
		virtual const std::vector<std::string> get_LongKey() const = 0;
		virtual const std::vector<std::string> get_ParamKey() const = 0;
		virtual const std::vector<std::string> get_GroupRateKey() const = 0;
		virtual const std::vector<std::string> get_AugmKey() const = 0;

		/////////////////////
		//////  MODEL  //////
		/////////////////////

		/// Whether the model uses contact rates.
		virtual bool has_contact_rates() const { return false; }

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

		/// Initialize advanced setup.
		virtual void init_advanced_setup() {}

		/// Initialize templated infecvitiy profile.
		void init_infec_profile_template(dict<std::string> x) { mcmc_->infec_profile.init_template(mcmc_->param, x); }

		/// Drop all augmented definitions of given individual.
		void clear_augm(Individual &member) { mcmc_->augm.clear_definitions(member); }

		/// User call. Initialize all augmented variables of given individual.
		/// If a variable has been set by the user, it's skipped.
		/// UPDATE: right now simulations assignments must be order independent.
		void init_missing_augm(Individual &member)
		{
			mcmc_->augm.init_missing_definitions(member);
			mcmc_->augm.user_reset(mcmc_->setup.prng(), member);
		}

		/// User call. Request discrete parameter definition.
		void request_discrete_param(size_t p) { mcmc_->param.request_discrete(p); }
		/// User call. Request discrete augmented definition.
		void request_discrete_augm(size_t a) { mcmc_->augm.request_discrete(a); }
		/// User call. Request discrete infectivity profile definition.
		void request_discrete_infec_profile() { mcmc_->infec_profile.request_discrete(); }

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		/// Sample a uniform-distributed variable between 0 and 1.
		double uniform() const { return mcmc_->setup.uniform(); }

		/// Sample a exponential-distributed variable between 0 and Inf.
		double exponential() const { return mcmc_->setup.exponential(); }

		/// Wheter hosuehold members i and j are equivalent in terms of rates.
		bool is_same_rate_partition(const Household &house, size_t i, size_t j) const
		{
			return mcmc_->is_same_rate_partition(house, i, j);
		}

		/// Get parameter value.
		double par(size_t p) const { return mcmc_->param.at(p).value(); }

		/// Return maximum number of observed infected household members.
		size_t get_max_hh_cases() const noexcept { return mcmc_->data_wide().get_max_hh_cases(); }
		/// Throw an exception if at least a member as missing status.
		void assert_no_missing_status() const { mcmc_->data_wide().assert_no_missing_status(); }

		double infec_profile_inst(const Individual &i, double t) const { return mcmc_->infec_profile.inst(&i, t); }
		double infec_profile_cumul(const Individual &i, double t) const { return mcmc_->infec_profile.cumul(&i, t); }

		/// Return reference to cached size_t value.
		size_t &cache_size(Household &house, size_t e) { return mcmc_->data_wide().cache_size(house, e); }
		/// Return cached size_t value.
		size_t cache_size(const Household &house, size_t e) const { return mcmc_->data_wide().cache_size(house, e); }

		/// Return cached size_t vector.
		std::vector<size_t> &cache_size_vector(Household &house, size_t e) { return mcmc_->data_wide().cache_size_vector(house, e); }
		const std::vector<size_t> &cache_size_vector(const Household &house, size_t e) const { return mcmc_->data_wide().cache_size_vector(house, e); }

		/// Return cached double vector.
		std::vector<double> &cache_double_vector(Household &house, size_t e) { return mcmc_->data_wide().cache_double_vector(house, e); }
		const std::vector<double> &cache_double_vector(const Household &house, size_t e) const { return mcmc_->data_wide().cache_double_vector(house, e); }

		/// Return cached double matrix.
		std::vector<std::vector<double>> &cache_double_matrix(Household &house, size_t e) { return mcmc_->data_wide().cache_double_matrix(house, e); }
		const std::vector<std::vector<double>> &cache_double_matrix(const Household &house, size_t e) const { return mcmc_->data_wide().cache_double_matrix(house, e); }

	private:
		MCMCBase *mcmc_;
	};
}
