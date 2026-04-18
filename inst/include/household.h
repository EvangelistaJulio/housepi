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
#include "house_storage.h"
#include "individual.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Household : class ***  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
	/// Household class.
	class Household
	{
		friend struct DataWide;
		friend struct ContactRate;
		friend class MCMCBase;
		friend class MCMC;

	public:
		//////////////////////
		//////  SETTER  //////
		//////////////////////

		/// Emplace a new individual in the household.
		void add_member(size_t index, const std::string &id, const Indiv::Wide *data_ptr);

		/// Set status of member 'i' to susceptible.
		void set_susceptible_status(size_t i) { helper_set_status_(i, status::SUSCEPTIBLE); }
		void set_susceptible_status(Individual &member) { set_susceptible_status(member.hhindex()); }

		/// Set status of member 'i' to infected.
		void set_infected_status(size_t i, int new_status = status::INFECTED)
		{
			if (!status::is_infected(new_status))
				raise_set_non_infec_status_(new_status);
			helper_set_status_(i, new_status);
		}
		void set_infected_status(Individual &member, int new_status = status::INFECTED)
		{
			set_infected_status(member.hhindex(), new_status);
		}

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		Individual &member(size_t i) { return members_[i]; }
		const Individual &member(size_t i) const { return members_[i]; }

		const std::vector<size_t> &infected_hhindexes() const noexcept { return infected_list_; }
		const std::vector<size_t> &susceptible_hhindexes() const noexcept { return healthy_list_; }
		const std::vector<size_t> &missing_hhindexes() const noexcept { return missing_list_; }

		size_t n_infected_status() const noexcept { return infected_list_.size(); }
		size_t n_susceptible_status() const noexcept { return healthy_list_.size(); }
		size_t n_missing_status() const noexcept { return missing_list_.size(); }

		/// Contact rate 'from' household member 'to member.
		/// No self-contact (set 0 for safety).
		double contact_rate(size_t from, size_t to) const
		{
			if (from == to)
				return 0.;
			if (!pcontact_rate_)
				return 1.;
			return pcontact_rate_->mtx[from][to].value();
		}
		/// Contact rate intra member partition (self-contact).
		double intra_partition_contact_rate(size_t i) const
		{
			if (!pcontact_rate_)
				return 1.;
			return pcontact_rate_->mtx[i][i].value();
		}

		double obs_log_lik() const { return obs_log_lik_.value(); }

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
		size_t size() const noexcept { return house_size_; }

	private:
		std::string house_id_ = "";						/// User house id
		size_t ridx_ = -1;								/// Restricted unique index
		std::vector<Individual> members_ = {};			/// Household members
		size_t house_size_ = 0;							/// Number of household members
		std::vector<size_t> infected_list_;				/// hh indexes of infected members
		std::vector<size_t> healthy_list_;				/// hh indexes of susceptible members
		std::vector<size_t> missing_list_;				/// hh indexes of members with missing status
		const House::Contact *pcontact_rate_ = nullptr; /// contact rates (see ContactRate)
		LogLike obs_log_lik_;							/// pointwise observation log-likelihood

		/// Change_status of individual i to new status (without controls).
		///	NOTE: this introduces augmented variables -> it requires RJMCMC in case of inference.
		void helper_set_status_(size_t i, int new_status);
		/// Raise error if status doesn't satisfy infected convention.
		[[noreturn]] void raise_set_non_infec_status_(int new_status) const;
	};
}
