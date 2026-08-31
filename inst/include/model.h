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

#define MAKE_HOUSEPI_CACHE_FORWARDER(NAME)                                                                              \
	auto &cache_##NAME(const Household &house, size_t e) noexcept { return mcmc_->data_wide().cache_##NAME(house, e); } \
	const auto &cache_##NAME(const Household &house, size_t e) const noexcept { return mcmc_->data_wide().cache_##NAME(house, e); }

	class ModelBase
	{
	public:
		/////////////////////
		///  CONSTRUCTOR  ///
		/////////////////////

		virtual ~ModelBase() = default;

		void bind_solver(MCMCBase *mcmc) noexcept { mcmc_ = mcmc; }

		virtual const std::vector<std::string> key_Data() const = 0;
		virtual const std::vector<std::string> key_HHmember() const = 0;
		virtual const std::vector<std::string> key_HHmembervec() const = 0;
		virtual const std::vector<std::string> key_HHdoublevec() const = 0;
		virtual const std::vector<std::string> key_HHdoublemtx() const = 0;
		virtual const std::vector<std::string> key_Long() const = 0;
		virtual const std::vector<std::string> key_Param() const = 0;
		virtual const std::vector<std::string> key_GroupRisk() const = 0;
		virtual const std::vector<std::string> key_Latent() const = 0;

		///////////////
		///  MODEL  ///
		///////////////

		/// Whether the model uses contact risks.
		virtual bool has_contact_risk() const { return false; }
		/// Whether the model requires an infectivity profile definition.
		virtual bool has_infec_profile() const { return false; }

		/// Recompute the household cached quantities consumed by the transmission model.
		virtual void update_house_transmission_cache(const Household & /*house*/) {}
		/// Same, restricted to a single individual.
		virtual void update_member_transmission_cache(const Household & /*house*/,
													  const Individual & /*indiv*/) {}
		/// Compute the transmission log-likelihood within a single household.
		virtual double compute_house_transm_log_lik(const Household & /*house*/) = 0;

		/// Inference specific initialization (e. g. cache variables, lookup tables, etc).
		virtual void init_model_for_evaluation() {}
		/// Inference specific initialization for each household.
		/// NOTE: it's run right after 'init_model_for_evaluation'.
		virtual void init_house_for_evaluation(const Household & /*house*/) {}

		/// Simulation specific initialization (e. g. cache variables, lookup tables, etc).
		virtual void init_model_for_simulation() {}
		/// Simulation for a single household.
		virtual void simulate_single_house(Household & /*house*/) = 0;

		/////////////////////////////////
		///  SETTER : MODEL REQUESTS  ///
		/////////////////////////////////

		/// Initialize model-specific definition requirements.
		virtual void init_model_requirements() {}

		/// Request that a baseline parameter and the relative risks scaling it share one reference.
		/// 'contacts' declares that contact risks scale the same term as the listed group risks.
		void request_risk_identifiability(size_t p, std::vector<size_t> risks, bool contacts = false)
		{
			mcmc_->param.request_risk_identifiability(p, std::move(risks), contacts);
		}

		/// Request discrete parameter definition.
		void request_param_discrete(std::vector<size_t> vec) { mcmc_->param.request_discrete(vec); }
		/// Request positive parameter definition.
		void request_param_positive(std::vector<size_t> vec) { mcmc_->param.request_positive(vec); }

		/// Request discrete latent definition.
		void request_latent_discrete(std::vector<size_t> vec) { mcmc_->augm.request_discrete(vec); }
		/// Request positive latent definition.
		void request_latent_positive(std::vector<size_t> vec) { mcmc_->augm.request_positive(vec); }
		/// Declare a latent key read by the transmission computation (included cache).
		void request_latent_transmission(std::vector<size_t> vec) { mcmc_->augm.request_transmission(vec); }

		/// Request discrete infectivity profile definition.
		void request_infec_profile_discrete() { mcmc_->infec_profile.request_discrete(); }

		/// Request discrete data definition.
		void request_data_discrete(std::vector<size_t> vec) { mcmc_->data_wide().request_discrete(vec); }
		/// Request data definition greater than zero.
		void request_data_positive(std::vector<size_t> vec) { mcmc_->data_wide().request_positive(vec); }

		/// Throw an exception if at least a member as missing status.
		void assert_no_missing_status() const { mcmc_->data_wide().assert_no_missing_status(); }

		//////////////////////////////
		///  SETTER : MODEL CALLS  ///
		//////////////////////////////

		/// Set an individual status to infected. Resolve status observation, if any.
		/// NOTE: call init_missing_latent to materialize status-specific definition block after.
		void set_infected(Household &house, Individual &indiv)
		{
			mcmc_->augm.simulate_infected_status(house, indiv, mcmc_->setup.prng());
		}
		/// Set an individual status to escaped status.
		void set_escaped(Household &house, Individual &indiv)
		{
			house.set_new_status_(indiv, status::ESCAPED);
		}

		/// Set numeric value of latent variable. If not allowed, raise exception.
		void set_latent_value(size_t a, const Individual &m, double x) { mcmc_->get_latent_varia(a, m).set(x); }

		/// Initialize all latent variables and observations.
		/// If a variable has been set by the user, it's skipped.
		/// UPDATE: right now simulations assignments must be order independent.
		void init_missing_latent(Individual &indiv)
		{
			mcmc_->augm.init_definitions_for_simulation(indiv);
			mcmc_->augm.draw_sample(indiv, mcmc_->setup.prng());
			mcmc_->augm.sample_obs(indiv, mcmc_->setup.prng());
		}

		////////////////
		///  GETTER  ///
		////////////////

		/// Sample a uniform-distributed variable between 0 and 1.
		double uniform() const { return mcmc_->setup.uniform(); }
		/// Sample a exponential-distributed variable between 0 and Inf.
		double exponential() const { return mcmc_->setup.exponential(); }

		/// Real-time integration step for simulation. Returns the user-configured.
		double real_time_step() const { return mcmc_->setup.time_step; }

		/// Get parameter value.
		double par(size_t p) const { return mcmc_->param.at(p).value(); }
		/// Get latent value for given individual.
		double aug(size_t a, const Individual &m) const { return mcmc_->aug(a, m); }
		/// Get wide-format data value for given individual.
		double data(size_t d, const Individual &m) const { return mcmc_->data(d, m); }
		/// Get group-risk value for given individual.
		double group_risk(size_t g, const Individual &m) const { return mcmc_->group_risk.value(g, m); }

		/// Names for error messages.
		std::string err_name(const Individual &m) const { return mcmc_->data_wide().err_name(m); }
		std::string err_name(const Household &h) const { return mcmc_->data_wide().err_name(h); }

		/// Get contact risk 'from' household member 'to' member.
		double contact_risk(const Household &house, MemberIndex from, MemberIndex to) const { return mcmc_->contact_risk.value(house, from, to); }
		/// Wheter hosuehold members i and j are equivalent in terms of risks.
		bool is_same_risk_partition(const Household &house, MemberIndex i, MemberIndex j) const
		{
			return mcmc_->is_same_risk_partition(house, i, j);
		}
		/// Get contact risk intra member partition (self-contact).
		double intra_partition_contact_risk(const Household &house, MemberIndex i) const { return mcmc_->contact_risk.intra_value(house, i); }

		/// Return maximum number of observed infected household members.
		size_t get_max_hh_cases() const noexcept { return mcmc_->data_wide().get_max_hh_cases(); }

		double infec_profile_inst(const Individual &m, double t) const { return mcmc_->infec_profile.inst(&m, t); }
		double infec_profile_cumul(const Individual &m, double t) const { return mcmc_->infec_profile.cumul(&m, t); }

		///////////////
		///  CACHE  ///
		///////////////

		/// Return cached member value.
		MAKE_HOUSEPI_CACHE_FORWARDER(member)
		/// Return cached member vector.
		MAKE_HOUSEPI_CACHE_FORWARDER(member_vector)
		/// Return cached double vector.
		MAKE_HOUSEPI_CACHE_FORWARDER(double_vector)
		/// Return cached double matrix.
		MAKE_HOUSEPI_CACHE_FORWARDER(double_matrix)

	private:
		MCMCBase *mcmc_;
	};
}

#undef MAKE_HOUSEPI_CACHE_FORWARDER
