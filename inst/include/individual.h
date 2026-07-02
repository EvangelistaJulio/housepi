/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>

#include "tools/messages.h"
#include "tools/status.h"

namespace housepi
{
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////// *** Individual : class *** /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	/// Class for individual household member.
	class Individual final
	{
		friend class Household; /// owns and mutates its members
		friend struct Index;	/// reads the opaque storage index_

	public:
		/////////////////////////
		//////  OPERATORS  //////
		/////////////////////////

		/// Equal to operator.
		bool operator==(const Individual &other) const noexcept { return index_ == other.index_; }
		/// Not equal to operator.
		bool operator!=(const Individual &other) const noexcept { return index_ != other.index_; }

		//////////////////////
		//////  GETTER  //////
		//////////////////////

		/// Unique identifier as provided in data.
		const std::string &id() const noexcept { return id_; }
		/// Associated status.
		int status() const noexcept { return status_; }
		/// Within household relative-index.
		/// UPDATE: this shall be dropped in the future
		size_t hhindex() const noexcept { return hhindex_; }

		/////////////////////
		//////  FLAGS  //////
		/////////////////////

		bool is_infected_status() const noexcept { return status::is_infected(status_); }
		bool is_escaped_status() const noexcept { return status::is_escaped(status_); }
		bool is_missing_status() const noexcept { return status::is_missing(status_); }

		////////////////////
		//////  SHOW  //////
		////////////////////

		std::string err_name() const { return "individual " + msg::val(id()); }

	private:
		size_t index_ = 0;	  /// housepi internal index
		size_t hhindex_ = 0;  /// household relative index
		std::string id_ = ""; /// user id

		/// status of the individual
		int status_ = status::MISSING;
	};

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// *** Index : struct *** ///////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	class Household;

	/// Index of Individual and Household. Client-Attorney idiom for lightweight viewers.
	struct Index
	{
		static size_t of(const Individual &member) noexcept { return member.index_; }
		static size_t of(const Household &house) noexcept; /// NOTE: defined in household.h
	};
}
