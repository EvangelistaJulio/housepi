/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

namespace housepi
{
	namespace status
	{
		constexpr int INFECTED = 1;	   /// Infected status (>=)
		constexpr int SUSCEPTIBLE = 0; /// Susceptible status (==)
		constexpr int MISSING = -1;	   /// Missing status (<=)

		/// Wheter individual status is infected.
		inline bool is_infected(int s) { return s >= INFECTED; }
		/// Wheter individual status is susceptible.
		inline bool is_susceptible(int s) { return s == SUSCEPTIBLE; }
		/// Wheter individual status is missing.
		inline bool is_missing(int s) { return s <= MISSING; }
	}
}
