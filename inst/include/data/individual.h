/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <cstdint>

#include "data/status.h"

//////////////////////////////  Individual  //////////////////////////////

namespace housepi
{
	using IndivIdx = std::uint32_t; /// Global individual index.
	using HouseIdx = std::uint32_t; /// Global household index.
	using MemberIdx = std::uint8_t; /// Member index within its household.

	/// Member of one household.
	///
	/// The user id lives in data::Wide, keyed by individual index, see Index::of().
	class Individual final
	{
		friend class Household; /// Owns and mutates its members.
		friend struct Index;	/// Reads the opaque storage indexes.

	public:
		//////  Operator  //////

		/// Whether they are different individuals.
		bool operator!=(const Individual &other) const noexcept { return indiv_index_ != other.indiv_index_; }

		//////  Getter  //////

		/// Status.
		StatusLvl status() const noexcept { return status_; }

		/// Whether status is infected.
		bool is_infected_status() const noexcept { return status::is_infected(status_); }
		/// Whether status is missing.
		bool is_missing_status() const noexcept { return status::is_missing(status_); }

	private:
		IndivIdx indiv_index_ = 0;			 /// Global individual index.
		MemberIdx member_index_ = 0;		 /// Index within its household.
		StatusLvl status_ = status::MISSING; /// Status.
	};

	//////////////////////////////  Index  //////////////////////////////

	class Household;

	/// Storage indexes of Individual and Household, a client-attorney for lightweight viewers.
	struct Index
	{
		/// Index of `indiv` within its household.
		static MemberIdx member(const Individual &indiv) noexcept { return indiv.member_index_; }
		/// Global index of `indiv`.
		static IndivIdx of(const Individual &indiv) noexcept { return indiv.indiv_index_; }
		/// Global index of `house`, defined in data/household.h.
		static HouseIdx of(const Household &house) noexcept;
	};
}
