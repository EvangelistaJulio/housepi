/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <limits>

#include "types.h"

namespace housepi::status
{
	constexpr StatusLevel ESCAPED = 0;	/// escaped status (==)
	constexpr StatusLevel INFECTED = 1; /// infected status (>= and != MISSING)
	/// missing status (==), the top level a status can hold
	constexpr StatusLevel MISSING = std::numeric_limits<StatusLevel>::max();
	constexpr StatusLevel MAX_INFECTED = MISSING - 1; /// highest infected status

	// NOTE: the predicates take a wide status so that they also read the status
	// column of an R data frame, where a missing status is NA_INTEGER (== INT_MIN):
	// it is out of the infected and escaped domains, hence a missing category.

	/// Wheter individual status is infected.
	inline bool is_infected(int s) { return s >= INFECTED && s != MISSING; }
	/// Wheter individual status is escaped.
	inline bool is_escaped(int s) { return s == ESCAPED; }
	/// Wheter individual status is missing.
	inline bool is_missing(int s) { return s == MISSING; }
	// Associate a category (0=infected, 1=escaped, 2=missing).
	inline size_t category(int s)
	{
		if (is_infected(s))
			return 0;
		if (is_escaped(s))
			return 1;
		return 2;
	};
}
