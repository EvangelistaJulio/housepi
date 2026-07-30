/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>

#include "tools/status.h"
#include "tools/types.h"

namespace housepi
{
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// *** Individual : class *** /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	/// Class for individual household member.
	/// NOTE: the user id lives in DataWide, keyed by the individual index (see Index::of).
	class alignas(8) Individual final
	{
		friend class Household; /// owns and mutates its members
		friend struct Index;	/// reads the opaque storage indexes

	public:
		/////////////////////////
		//////  OPERATORS  //////
		/////////////////////////

		/// Equal to operator.
		bool operator==(const Individual &other) const noexcept { return indiv_index_ == other.indiv_index_; }
		/// Not equal to operator.
		bool operator!=(const Individual &other) const noexcept { return indiv_index_ != other.indiv_index_; }

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		/// Associated status.
		int status() const noexcept { return status_; }

		/////////////////////
		//////  FLAGS  //////
		/////////////////////

		bool is_infected_status() const noexcept { return status::is_infected(status_); }
		bool is_escaped_status() const noexcept { return status::is_escaped(status_); }
		bool is_missing_status() const noexcept { return status::is_missing(status_); }

	private:
		IndivIndex indiv_index_ = 0;   /// housepi internal individual index
		MemberIndex member_index_ = 0; /// index within its own household

		/// status of the individual
		int status_ = status::MISSING;
	};

	static_assert(sizeof(Individual) % 8 == 0, "Individual must be 8-byte padded.");

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// *** Index : struct *** ///////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	class Household;

	/// Index of Individual and Household. Client-Attorney idiom for lightweight viewers.
	struct Index
	{
		/// Global individual index.
		static size_t of(const Individual &member) noexcept { return member.indiv_index_; }
		/// Global household index.
		static size_t of(const Household &house) noexcept; /// NOTE: defined in household.h
		/// Index of an individual within its own household.
		static size_t member(const Individual &member) noexcept { return member.member_index_; }
	};
}
