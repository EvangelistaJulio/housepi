/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
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

		/// Auxiliary iterator of household member indexes with same status.
		struct IndexRange_
		{
			const size_t *data_;
			size_t size_;
			const size_t *begin() const noexcept { return data_; }
			const size_t *end() const noexcept { return data_ + size_; }
			size_t size() const noexcept { return size_; }
			bool empty() const noexcept { return size_ == 0; }
			const size_t &operator[](size_t i) const noexcept { return data_[i]; }
			const size_t &front() const noexcept { return data_[0]; }
		};

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
		void set_infected_status(size_t i, int new_status = status::INFECTED) { helper_set_status_(i, new_status); }
		void set_infected_status(Individual &member, int new_status = status::INFECTED) { set_infected_status(member.hhindex(), new_status); }

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		Individual &member(size_t i) noexcept { return members_[i]; }
		const Individual &member(size_t i) const noexcept { return members_[i]; }

		/// Iterator for household indexes of members with infected status.
		IndexRange_ infected_hhindexes() const noexcept { return {status_vec_.data(), ninf_}; }
		/// Iterator for household indexes of members with susceptible status.
		IndexRange_ susceptible_hhindexes() const noexcept { return {status_vec_.data() + ninf_, nsus_}; }
		/// Iterator for household indexes of members with missing status.
		IndexRange_ missing_hhindexes() const noexcept { return {status_vec_.data() + ninf_ + nsus_, house_size_ - ninf_ - nsus_}; }

		/// Number of members with infected status.
		size_t n_infected_status() const noexcept { return ninf_; }
		/// Number of members with susceptible status.
		size_t n_susceptible_status() const noexcept { return nsus_; }
		/// Number of members with missing status.
		size_t n_missing_status() const noexcept { return house_size_ - ninf_ - nsus_; }

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

		/// Value of household transmission log-likelihood.
		double transm_log_lik() const { return transm_log_lik_.value(); }

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

		/// Number of household members.
		size_t size() const noexcept { return house_size_; }
		/// Index of equivalence class.
		size_t eq_class() const noexcept { return eq_class_; }
		/// Household ID as in data.
		const std::string &house_id() const noexcept { return house_id_; }
		/// Name for error messages.
		std::string err_name() const { return "household '" + house_id_ + "'"; }

	private:
		std::string house_id_ = "";						/// user house id
		size_t eq_class_ = -1;							/// equivalence class index
		std::vector<Individual> members_ = {};			/// household members
		std::vector<size_t> status_vec_ = {};			/// partitioned: [infected | susceptible | missing]
		size_t house_size_ = 0;							/// number of household members
		size_t ninf_ = 0;								/// number of infected members
		size_t nsus_ = 0;								/// number of susceptible members
		const House::Contact *pcontact_rate_ = nullptr; /// contact rates (see ContactRate)
		LogDensity transm_log_lik_;						/// transmission log-likelihood

		/// Change_status of individual i to new status (without controls).
		/// NOTE: Using 'partition refinement' to track member statuses.
		void helper_set_status_(size_t i, int new_status);
	};
}
