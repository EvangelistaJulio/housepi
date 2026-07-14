/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
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

		virtual ~ModelBase() = default;

		void bind_solver(MCMCBase *mcmc) noexcept { mcmc_ = mcmc; }

		virtual const std::vector<std::string> get_DataKey() const = 0;
		virtual const std::vector<std::string> get_HHuintKey() const = 0;
		virtual const std::vector<std::string> get_HHuintvecKey() const = 0;
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
		virtual void update_house_cache_dependecies(Household & /*house*/) {}
		/// Update cached augmeted dependecies for single household member.
		virtual void update_member_cache_dependecies(Household & /*house*/, Individual & /*member*/) {}
		/// Compute the transmission log-likelihood within a single household.
		virtual double compute_house_transm_log_lik(Household & /*house*/) = 0;

		/// Inference specific initialization (e. g. cache variables, lookup tables, etc).
		virtual void init_model_for_evaluation() {}
		/// Inference specific initialization for each household.
		/// NOTE: it's run right after 'init_model_for_evaluation'.
		virtual void init_house_for_evaluation(Household & /*house*/) {}

		/// Simulation specific initialization (e. g. cache variables, lookup tables, etc).
		virtual void init_model_for_simulation() {}
		/// Simulation for a single household.
		virtual void simulate_single_house(Household & /*house*/) = 0;

		//////////////////////
		//////  SETTER  //////
		//////////////////////

		/// Initialize advanced setup.
		virtual void init_advanced_setup() {}

		/// Initialize templated infecvitiy profile.
		void init_infec_profile_template(dict<std::string> x) { mcmc_->infec_profile.init_template(mcmc_->translator, x); }

		/// Throw an exception if at least a member as missing status.
		void assert_no_missing_status() const { mcmc_->data_wide().assert_no_missing_status(); }

		//////////////////////////////////////
		//////  SETTER : MODELER CALLS  //////
		//////////////////////////////////////

		/// Request discrete parameter definition.
		void request_param_discrete(size_t p) { mcmc_->param.request_discrete(p); }
		/// Request discrete augmented definition.
		void request_augm_discrete(size_t a) { mcmc_->augm.request_discrete(a); }
		/// Request discrete infectivity profile definition.
		void request_infec_profile_discrete() { mcmc_->infec_profile.request_discrete(); }

		/// Request discrete data definition.
		void request_data_discrete(size_t d) { mcmc_->data_wide().request_discrete(d); }
		/// Request data definition greater than zero.
		void request_data_positive(size_t d) { mcmc_->data_wide().request_positive(d); }

		/// Set a member status to infected. Resolve status observation, if any.
		/// NOTE: call init_missing_augm to materialize status-specific definition block after.
		void set_infected(Household &house, Individual &member)
		{
			mcmc_->augm.simulate_infected_status(house, member, mcmc_->setup.prng());
		}
		/// Set a member status to escaped status.
		void set_escaped(Household &house, Individual &member)
		{
			house.set_new_status_(member, status::ESCAPED);
		}

		/// Drop all augmented definitions of given individual.
		void clear_augm(Individual &member) { mcmc_->augm.clear_definitions(member); }

		/// Set numeric value of augmented variable. If not allowed, raise exception.
		void set_augm_value(size_t a, const Individual &m, double x) { mcmc_->get_augm_varia(a, m).set(x); }

		/// Initialize all augmented variables and observations.
		/// If a variable has been set by the user, it's skipped.
		/// UPDATE: right now simulations assignments must be order independent.
		void init_missing_augm(Individual &member)
		{
			mcmc_->augm.init_definitions_for_simulation(member);
			mcmc_->augm.restart(member, mcmc_->setup.prng());
			mcmc_->augm.sample_obs(member, mcmc_->setup.prng());
		}

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		/// Sample a uniform-distributed variable between 0 and 1.
		double uniform() const { return mcmc_->setup.uniform(); }
		/// Sample a exponential-distributed variable between 0 and Inf.
		double exponential() const { return mcmc_->setup.exponential(); }

		/// Real-time integration step for simulation. Returns the user-configured.
		double real_time_step() const { return mcmc_->setup.time_step; }

		/// Get parameter value.
		double par(size_t p) const { return mcmc_->param.at(p).value(); }
		/// Get augmented value for given individual.
		double aug(size_t a, const Individual &m) const { return mcmc_->aug(a, m); }
		/// Get wide-format data value for given individual.
		double data(size_t d, const Individual &m) const { return mcmc_->data(d, m); }
		/// Get group-rate value for given individual.
		double group_rate(size_t g, const Individual &m) const { return mcmc_->group_rate.value(g, m); }

		/// Get contact rate 'from' household member 'to' member.
		double contact_rate(const Household &house, size_t from, size_t to) const { return mcmc_->contact_rate.value(house, from, to); }
		/// Wheter hosuehold members i and j are equivalent in terms of rates.
		bool is_same_rate_partition(const Household &house, size_t i, size_t j) const
		{
			return mcmc_->is_same_rate_partition(house, i, j);
		}
		/// Get contact rate intra member partition (self-contact).
		double intra_partition_contact_rate(const Household &house, size_t i) const { return mcmc_->contact_rate.intra_value(house, i); }

		/// Return maximum number of observed infected household members.
		size_t get_max_hh_cases() const noexcept { return mcmc_->data_wide().get_max_hh_cases(); }

		double infec_profile_inst(const Individual &m, double t) const { return mcmc_->infec_profile.inst(&m, t); }
		double infec_profile_cumul(const Individual &m, double t) const { return mcmc_->infec_profile.cumul(&m, t); }

		/// Return reference to cached unsigned value.
		unsigned &cache_uint(Household &house, size_t e) noexcept { return mcmc_->data_wide().cache_uint(house, e); }
		/// Return cached unsigned value.
		unsigned cache_uint(const Household &house, size_t e) const noexcept { return mcmc_->data_wide().cache_uint(house, e); }

		/// Return cached unsigned vector.
		std::vector<unsigned> &cache_uint_vector(Household &house, size_t e) noexcept { return mcmc_->data_wide().cache_uint_vector(house, e); }
		const std::vector<unsigned> &cache_uint_vector(const Household &house, size_t e) const noexcept { return mcmc_->data_wide().cache_uint_vector(house, e); }

		/// Return cached double vector.
		std::vector<double> &cache_double_vector(Household &house, size_t e) noexcept { return mcmc_->data_wide().cache_double_vector(house, e); }
		const std::vector<double> &cache_double_vector(const Household &house, size_t e) const noexcept { return mcmc_->data_wide().cache_double_vector(house, e); }

		/// Return cached double matrix.
		std::vector<std::vector<double>> &cache_double_matrix(Household &house, size_t e) noexcept { return mcmc_->data_wide().cache_double_matrix(house, e); }
		const std::vector<std::vector<double>> &cache_double_matrix(const Household &house, size_t e) const noexcept { return mcmc_->data_wide().cache_double_matrix(house, e); }

	private:
		MCMCBase *mcmc_;
	};
}
