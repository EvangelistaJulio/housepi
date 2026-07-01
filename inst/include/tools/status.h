/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

namespace housepi::status
{
	constexpr int INFECTED = 1; /// infected status (>=)
	constexpr int ESCAPED = 0;	/// escaped status (==)
	constexpr int MISSING = -1; /// missing status (<=)

	/// Wheter individual status is infected.
	inline bool is_infected(int s) { return s >= INFECTED; }
	/// Wheter individual status is escaped.
	inline bool is_escaped(int s) { return s == ESCAPED; }
	/// Wheter individual status is missing.
	inline bool is_missing(int s) { return s <= MISSING; }
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
