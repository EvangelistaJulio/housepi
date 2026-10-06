/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>
#include <string>

#include "data/individual.h"
#include "mcmc/base.h"

//////////////////////////////  ModelBase  //////////////////////////////

namespace housepi
{
	class Individual;
	class Household;

#define MAKE_HOUSEPI_CACHE_FORWARDER(NAME)                                    \
	auto &cache_##NAME(const Household &house, size_t e) noexcept             \
	{                                                                         \
		return mcmc_->data.state.cache_##NAME(house, e);                      \
	}                                                                         \
	const auto &cache_##NAME(const Household &house, size_t e) const noexcept \
	{                                                                         \
		return mcmc_->data.state.cache_##NAME(house, e);                      \
	}

	/// Transmission model: per-household log-likelihood and simulation hooks.
	///
	/// - built-in and custom models derive from it through the `HOUSEPI_*` macros;
	/// - the accessors below read the engine the model is bound to.
	class ModelBase
	{
	public:
		//////  Constructor  //////

		virtual ~ModelBase() = default;

		/// Bind the engine `mcmc` the model reads.
		void bind_solver(MCMCBase *mcmc) noexcept { mcmc_ = mcmc; }

		//////  Model-declared keys  //////

		virtual const std::vector<std::string> key_Data() const = 0;
		virtual const std::vector<std::string> key_HHmember() const = 0;
		virtual const std::vector<std::string> key_HHmembervec() const = 0;
		virtual const std::vector<std::string> key_HHdoublevec() const = 0;
		virtual const std::vector<std::string> key_HHdoublemtx() const = 0;
		virtual const std::vector<std::string> key_Long() const = 0;
		virtual const std::vector<std::string> key_Param() const = 0;
		virtual const std::vector<std::string> key_GroupRisk() const = 0;
		virtual const std::vector<std::string> key_Latent() const = 0;

		//////  Model  //////

		/// Whether the model uses contact risks.
		virtual bool has_contact_risk() const { return false; }
		/// Whether the model requires an infectivity profile definition.
		virtual bool has_infec_profile() const { return false; }

		/// Recompute the cached quantities of `house` read by the transmission model.
		virtual void update_house_transmission_cache(const Household & /*house*/) {}
		/// Recompute the cached quantities of member `indiv` read by the transmission model.
		virtual void update_member_transmission_cache(const Household & /*house*/,
													  const Individual & /*indiv*/) {}
		/// Compute the transmission log-likelihood of `house`.
		virtual double compute_house_transm_log_lik(const Household & /*house*/) = 0;

		/// Initialize the model for evaluation (cache variables, lookup tables).
		virtual void init_model_for_evaluation() {}
		/// Initialize household `house` for evaluation.
		/// - runs right after init_model_for_evaluation().
		virtual void init_house_for_evaluation(const Household & /*house*/) {}

		/// Initialize the model for simulation (cache variables, lookup tables).
		virtual void init_model_for_simulation() {}
		/// Simulate household `house`.
		virtual void simulate_single_house(Household & /*house*/) = 0;

		//////  Setter: model requests  //////

		/// Declare the model definition requirements.
		virtual void init_model_requirements() {}

		/// Request that baseline parameter `p` and group risks `risks` share one reference.
		/// @param contacts Whether contact risks scale the same term as `risks`.
		void request_risk_identifiability(size_t p, std::vector<size_t> risks, bool contacts = false)
		{
			mcmc_->risks.request_identifiability(p, std::move(risks), contacts);
		}

		/// Declare the group risks `vec` a rule may stratify by infected status.
		/// @note Leave out risks read before a member's status is drawn in simulation.
		void request_risk_status_stratified(std::vector<size_t> vec)
		{
			mcmc_->risks.request_status_stratified(std::move(vec));
		}

		/// Request discrete parameters `vec`.
		void request_param_discrete(std::vector<size_t> vec) { mcmc_->param.request_discrete(vec); }
		/// Request positive parameters `vec`.
		void request_param_positive(std::vector<size_t> vec) { mcmc_->param.request_positive(vec); }

		/// Request discrete latents `vec`.
		void request_latent_discrete(std::vector<size_t> vec) { mcmc_->augm.request_discrete(vec); }
		/// Request positive latents `vec`.
		void request_latent_positive(std::vector<size_t> vec) { mcmc_->augm.request_positive(vec); }
		/// Declare latents `vec` as generated and read by transmission, cache included.
		/// A latent listed here takes no model: transmission gives its density.
		void request_latent_transmission(std::vector<size_t> vec) { mcmc_->augm.request_transmission(vec); }

		/// Request a discrete infectivity profile.
		void request_infec_profile_discrete() { mcmc_->infec_profile.request_discrete(); }

		/// Request discrete data columns `vec`.
		void request_data_discrete(std::vector<size_t> vec) { mcmc_->data.wide().request_discrete(vec); }
		/// Request positive data columns `vec`.
		void request_data_positive(std::vector<size_t> vec) { mcmc_->data.wide().request_positive(vec); }
		/// Request data columns `vec` with a value for every member, escaped ones included.
		void request_data_every_member(std::vector<size_t> vec) { mcmc_->data.wide().request_everyone(vec); }

		/// Throw if some member has missing status.
		void assert_no_missing_status() const { mcmc_->data.wide().assert_no_missing_status(); }

		//////  Setter: model calls  //////

		/// Set member `indiv` infected, drawing its status observation, if any.
		/// @note Call init_missing_latent() next, to build the definitions of its status.
		void set_infected(Household &house, Individual &indiv)
		{
			mcmc_->augm.simulate_infected_status(house, indiv, mcmc_->setup.prng());
		}
		/// Set member `indiv` escaped.
		void set_escaped(Household &house, Individual &indiv)
		{
			house.set_new_status_(indiv, status::ESCAPED);
		}

		/// Set latent `a` of member `m` to `x`; throw if not allowed.
		void set_latent_value(size_t a, const Individual &m, double x) { mcmc_->augm.at(a, m).set(x); }

		/// Initialize the latents and observations of member `indiv`.
		/// - latent set by the user: skipped.
		/// @todo Allow order-dependent assignments in simulation.
		void init_missing_latent(Individual &indiv)
		{
			mcmc_->augm.init_definitions_for_simulation(indiv);
			mcmc_->augm.draw_sample_for_simulation(indiv, mcmc_->setup.prng());
			mcmc_->augm.sample_obs(indiv, mcmc_->setup.prng());
		}

		//////  Getter  //////

		/// Draw a uniform variate on [0, 1].
		double uniform() const { return mcmc_->setup.uniform(); }
		/// Draw a unit-rate exponential variate.
		double exponential() const { return mcmc_->setup.exponential(); }

		/// User-configured time step of the simulation.
		double real_time_step() const { return mcmc_->setup.time_step; }

		/// Value of parameter `p`.
		double par(size_t p) const { return mcmc_->param.at(p).value(); }
		/// Latent `a` of member `m`.
		double aug(size_t a, const Individual &m) const { return mcmc_->aug(a, m); }
		/// Datum `d` of member `m`.
		double data(size_t d, const Individual &m) const { return mcmc_->data.at(d, m).value(); }
		/// Group risk `g` of member `m`.
		double group_risk(size_t g, const Individual &m) const { return mcmc_->risks.group().value(g, m); }

		/// Name of member `m`, for error messages.
		std::string str(const Individual &m) const { return mcmc_->data.wide().str(m); }
		/// Name of household `h`, for error messages.
		std::string str(const Household &h) const { return mcmc_->data.wide().str(h); }

		/// Contact risk from member `from` to member `to` of `house`.
		double contact_risk(const Household &house, MemberIdx from, MemberIdx to) const
		{
			return mcmc_->risks.contact().value(house, from, to);
		}
		/// Whether members `i` and `j` of `house` are equivalent according to risks.
		bool is_same_risk_partition(const Household &house, MemberIdx i, MemberIdx j) const
		{
			return mcmc_->risks.same_partition(house, i, j);
		}
		/// Contact risk within the partition of member `i` (self-contact).
		double intra_partition_contact_risk(const Household &house, MemberIdx i) const
		{
			return mcmc_->risks.contact().intra_value(house, i);
		}

		/// Maximum number of infected members observed in one household.
		size_t get_max_hh_cases() const noexcept { return mcmc_->data.wide().get_max_hh_cases(); }

		/// Instantaneous infectivity profile of member `m` at time `t`.
		double infec_profile_inst(const Individual &m, double t) const
		{
			return mcmc_->infec_profile.inst(&m, t);
		}
		/// Cumulative infectivity profile of member `m` at time `t`.
		double infec_profile_cumul(const Individual &m, double t) const
		{
			return mcmc_->infec_profile.cumul(&m, t);
		}

		//////  Cache  //////

		/// Cached member value.
		MAKE_HOUSEPI_CACHE_FORWARDER(member)
		/// Cached member vector.
		MAKE_HOUSEPI_CACHE_FORWARDER(member_vector)
		/// Cached double vector.
		MAKE_HOUSEPI_CACHE_FORWARDER(double_vector)
		/// Cached double matrix.
		MAKE_HOUSEPI_CACHE_FORWARDER(double_matrix)

	private:
		MCMCBase *mcmc_; /// Engine the model reads.
	};
}

#undef MAKE_HOUSEPI_CACHE_FORWARDER
