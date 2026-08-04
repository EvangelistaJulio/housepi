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

		/// Auxiliary iterator for individuals with same status.
		struct ConstIndivRange_
		{
			const Household *house_;
			const MemberIndex *data_;
			size_t size_;
			struct Iterator_
			{
				const Household *house_;
				const MemberIndex *ptr_;
				const Individual &operator*() const noexcept { return house_->indivs_[*ptr_]; }
				Iterator_ &operator++() noexcept
				{
					++ptr_;
					return *this;
				}
				bool operator!=(const Iterator_ &other) const noexcept { return ptr_ != other.ptr_; }
			};
			Iterator_ begin() const noexcept { return {house_, data_}; }
			Iterator_ end() const noexcept { return {house_, data_ + size_}; }
			size_t size() const noexcept { return size_; }
			bool empty() const noexcept { return size_ == 0; }
			const Individual &operator[](size_t k) const noexcept { return house_->indivs_[data_[k]]; }
		};

	public:
		//////////////////////
		//////  SETTER  //////
		//////////////////////

		/// Emplace a new individual in the household.
		void add_indiv(IndivIndex index, StatusLevel status);

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		/// Individual with member index i.
		Individual &at(MemberIndex i) noexcept { return indivs_[i]; }
		const Individual &at(MemberIndex i) const noexcept { return indivs_[i]; }

		/// Iterator for member indexes of members with infected status.
		MemberRange_ infected_members() const noexcept { return {status_vec_.data(), ninf_}; }
		/// Iterator for member indexes of members with escaped status.
		MemberRange_ escaped_members() const noexcept { return {status_vec_.data() + ninf_, nsus_}; }
		/// Iterator for member indexes of members with missing status.
		MemberRange_ missing_members() const noexcept { return {status_vec_.data() + ninf_ + nsus_, n_missing_status()}; }

		/// Iterator for individuals with infected status.
		ConstIndivRange_ infected_indivs() const noexcept { return {this, status_vec_.data(), ninf_}; }
		/// Iterator for individuals with escaped status.
		ConstIndivRange_ escaped_indivs() const noexcept { return {this, status_vec_.data() + ninf_, nsus_}; }
		/// Iterator for individuals with missing status.
		ConstIndivRange_ missing_indivs() const noexcept { return {this, status_vec_.data() + ninf_ + nsus_, n_missing_status()}; }

		/// Number of members with infected status.
		MemberIndex n_infected_status() const noexcept { return ninf_; }
		/// Number of members with escaped status.
		MemberIndex n_escaped_status() const noexcept { return nsus_; }
		/// Number of members with missing status.
		MemberIndex n_missing_status() const noexcept { return hhsize_ - ninf_ - nsus_; }

		/////////////////////////
		//////  ITERATORS  //////
		/////////////////////////

		auto begin() noexcept { return indivs_.begin(); }
		auto end() noexcept { return indivs_.end(); }
		auto begin() const noexcept { return indivs_.begin(); }
		auto end() const noexcept { return indivs_.end(); }

		////////////////////
		//////  SHOW  //////
		////////////////////

		/// Number of household members.
		MemberIndex size() const noexcept { return hhsize_; }
		/// Index of equivalence class.
		HouseIndex eq_class() const noexcept { return eq_class_; }

	private:
		HouseIndex house_index_ = 0;			   /// housepi internal household index
		HouseIndex eq_class_ = na::UINDEX;		   /// equivalence class index
		MemberIndex hhsize_ = 0;				   /// number of household members
		MemberIndex ninf_ = 0;					   /// number of infected members
		MemberIndex nsus_ = 0;					   /// number of escaped members
		std::vector<Individual> indivs_ = {};	   /// household members
		std::vector<MemberIndex> status_vec_ = {}; /// partitioned: [infected | escaped | missing]

		/// Change status of individual i to new status.
		/// NOTE: Using 'partition refinement' to track member statuses.
		void set_new_status_(Individual &indiv, StatusLevel new_status);
	};

	/// Index of Household (see individual.h).
	inline HouseIndex Index::of(const Household &house) noexcept { return house.house_index_; }
}
