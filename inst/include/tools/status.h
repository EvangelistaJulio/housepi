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
	/// escaped status (==)
	constexpr StatusLvl ESCAPED = 0;
	/// infected status (>= and != MISSING)
	constexpr StatusLvl INFECTED = 1;
	/// missing status (==), the top level a status can hold
	constexpr StatusLvl MISSING = std::numeric_limits<StatusLvl>::max();
	/// largest infected status
	constexpr StatusLvl MAX_INFECTED = MISSING - 1;

	/// Whether individual status is infected.
	inline bool is_infected(int s) { return s >= INFECTED && s != MISSING; }
	/// Whether individual status is escaped.
	inline bool is_escaped(int s) { return s == ESCAPED; }
	/// Whether individual status is missing.
	inline bool is_missing(int s) { return s == MISSING; }
	/// Associate a category (0=infected, 1=escaped, 2=missing).
	inline size_t category(int s) { return is_infected(s) ? 0 : (is_escaped(s) ? 1 : 2); };
}
