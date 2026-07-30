/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>

#include "tools/nan.h"
#include "tools/types.h"
#include "individual.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Household : class ***  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
	/// Household class.
	/// NOTE: the user house id lives in DataWide, keyed by the household index (see Index::of).
	class Household
	{
		friend struct Index;	 /// reads the opaque storage index
		friend struct DataWide;	 /// owns and mutates its households
		friend class ModelBase;	 /// changes member status
		friend struct Augmented; /// changes member status

		/// Auxiliary iterator for member indexes of members with same status.
		struct MemberRange_
		{
			const MemberIndex *data_;
			size_t size_;
			const MemberIndex *begin() const noexcept { return data_; }
			const MemberIndex *end() const noexcept { return data_ + size_; }
			size_t size() const noexcept { return size_; }
			bool empty() const noexcept { return size_ == 0; }
			const MemberIndex &operator[](size_t i) const noexcept { return data_[i]; }
			const MemberIndex &front() const noexcept { return data_[0]; }
		};

	public:
		//////////////////////
		//////  SETTER  //////
		//////////////////////

		/// Emplace a new individual in the household.
		void add_member(IndivIndex index, int status);

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		Individual &member(size_t i) noexcept { return members_[i]; }
		const Individual &member(size_t i) const noexcept { return members_[i]; }

		/// Iterator for member indexes of members with infected status.
		MemberRange_ infected_members() const noexcept { return {status_vec_.data(), ninf_}; }
		/// Iterator for member indexes of members with escaped status.
		MemberRange_ escaped_members() const noexcept { return {status_vec_.data() + ninf_, nsus_}; }
		/// Iterator for member indexes of members with missing status.
		MemberRange_ missing_members() const noexcept { return {status_vec_.data() + ninf_ + nsus_, n_missing_status()}; }

		/// Number of members with infected status.
		size_t n_infected_status() const noexcept { return ninf_; }
		/// Number of members with escaped status.
		size_t n_escaped_status() const noexcept { return nsus_; }
		/// Number of members with missing status.
		size_t n_missing_status() const noexcept { return static_cast<size_t>(hhsize_ - ninf_ - nsus_); }

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
		size_t size() const noexcept { return hhsize_; }
		/// Index of equivalence class.
		size_t eq_class() const noexcept { return eq_class_; }

	private:
		HouseIndex house_index_ = 0;			   /// housepi internal household index
		HouseIndex eq_class_ = na::UINDEX;		   /// equivalence class index
		MemberIndex hhsize_ = 0;				   /// number of household members
		MemberIndex ninf_ = 0;					   /// number of infected members
		MemberIndex nsus_ = 0;					   /// number of escaped members
		std::vector<Individual> members_ = {};	   /// household members
		std::vector<MemberIndex> status_vec_ = {}; /// partitioned: [infected | escaped | missing]

		/// Change status of individual i to new status.
		/// NOTE: Using 'partition refinement' to track member statuses.
		void set_new_status_(Individual &member, int new_status);
	};

	static_assert(sizeof(Household) % 8 == 0, "Household must be 8-byte padded.");

	/// Index of Household (see individual.h).
	inline size_t Index::of(const Household &house) noexcept { return house.house_index_; }
}
