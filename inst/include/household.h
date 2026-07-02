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

#include "tools/messages.h"
#include "tools/nan.h"
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
		friend struct Index;	 /// reads the opaque storage index_
		friend struct DataWide;	 /// owns and mutates its households
		friend class ModelBase;	 /// changes member status
		friend struct Augmented; /// changes member status

		/// Auxiliary iterator for household indexes of member with same status.
		struct HouseStatusRange_
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
		void add_member(size_t index, const std::string &id, int status);

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		Individual &member(size_t i) noexcept { return members_[i]; }
		const Individual &member(size_t i) const noexcept { return members_[i]; }

		/// Iterator for household indexes of members with infected status.
		HouseStatusRange_ infected_hhindexes() const noexcept { return {status_vec_.data(), ninf_}; }
		/// Iterator for household indexes of members with escaped status.
		HouseStatusRange_ escaped_hhindexes() const noexcept { return {status_vec_.data() + ninf_, nsus_}; }
		/// Iterator for household indexes of members with missing status.
		HouseStatusRange_ missing_hhindexes() const noexcept { return {status_vec_.data() + ninf_ + nsus_, house_size_ - ninf_ - nsus_}; }

		/// Number of members with infected status.
		size_t n_infected_status() const noexcept { return ninf_; }
		/// Number of members with escaped status.
		size_t n_escaped_status() const noexcept { return nsus_; }
		/// Number of members with missing status.
		size_t n_missing_status() const noexcept { return house_size_ - ninf_ - nsus_; }

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
		std::string err_name() const { return "household " + msg::val(house_id_); }

	private:
		std::string house_id_ = "";			   /// user house id
		size_t index_ = 0;					   /// housepi internal household index
		size_t eq_class_ = na::INDEX;		   /// equivalence class index
		std::vector<Individual> members_ = {}; /// household members
		std::vector<size_t> status_vec_ = {};  /// partitioned: [infected | escaped | missing]
		size_t house_size_ = 0;				   /// number of household members
		size_t ninf_ = 0;					   /// number of infected members
		size_t nsus_ = 0;					   /// number of escaped members

		/// Change status of individual i to new status.
		/// NOTE: Using 'partition refinement' to track member statuses.
		void set_new_status_(Individual &member, int new_status);
	};

	/// Index of Household (see individual.h).
	inline size_t Index::of(const Household &house) noexcept { return house.index_; }
}
