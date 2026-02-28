/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "function.h"
#include "individual.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Household : class ***  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
{
	/// Household class.
	class Household
	{
		friend struct DataWide;
		friend struct ContactRate;
		friend class MCMCBase;
		friend class MCMC;

	public:
		///////////////////////////
		//////  CONSTRUCTOR  //////
		///////////////////////////

		Household() = default;
		~Household() = default;

		//////////////////////
		//////  SETTER  //////
		//////////////////////

		/// Emplace a new individual in the household.
		void init_member(size_t index, const std::string &id,
						 int status, const std::vector<double> &data);

		void set_susceptible_status(Individual &member) { helper_set_status_(member, status::SUSCEPTIBLE); }
		void set_infected_status(Individual &member) { helper_set_status_(member, status::INFECTED); }
		void set_infected_status(Individual &member, int new_status)
		{
			if (!status::is_infected(new_status))
				throw std::invalid_argument(
					"Tried to set as infected status " + std::to_string(new_status) +
					" to " + err_name() + ".\nInfected status must be an integer greater than 0.");
			helper_set_status_(member, new_status);
		}

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		Individual &member(size_t i) { return members_.at(i); }
		const Individual &member(size_t i) const { return members_.at(i); }

		const std::vector<size_t> &infected_hhindexes() const noexcept { return infected_list_; }
		const std::vector<size_t> &susceptible_hhindexes() const noexcept { return healthy_list_; }
		const std::vector<size_t> &missing_hhindexes() const noexcept { return missing_list_; }

		size_t n_infected_status() const noexcept { return infected_list_.size(); }
		size_t n_susceptible_status() const noexcept { return healthy_list_.size(); }
		size_t n_missing_status() const noexcept { return missing_list_.size(); }

		double contact_rate(size_t from, size_t to) const { return contact_rates_[from][to].value(); }

		/////////////////////////
		//////  ITERATORS  //////
		/////////////////////////

		auto begin() noexcept { return members_.begin(); }
		auto end() noexcept { return members_.end(); }
		auto begin() const noexcept { return members_.begin(); }
		auto end() const noexcept { return members_.end(); }

		////////////////////
		//////  SHOW  //////
		////////////////////

		const std::string &house_id() const noexcept { return house_id_; }
		std::string err_name() const { return "household '" + house_id_ + "'"; }
		size_t size() const noexcept { return members_.size(); }

	private:
		std::string house_id_ = "";							/// User house id
		size_t restricted_ = -1;							/// Restricted unique index
		std::vector<Individual> members_ = {};				/// Household members
		std::vector<size_t> infected_list_;					/// hh indexes of infected members
		std::vector<size_t> healthy_list_;					/// hh indexes of susceptible members
		std::vector<size_t> missing_list_;					/// hh indexes of members with missing status
		std::vector<std::vector<Attribute>> contact_rates_; /// contact rates
		LogLike obs_log_lik_;								/// pointwise observation log-likelihood

		/// Change_status of individual i to new status (without controls).
		///	NOTE: this introduces augmented variables -> it requires RJMCMC in case of inference.
		void helper_set_status_(Individual &member, int new_status);
	};
}
