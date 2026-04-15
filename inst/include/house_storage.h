/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>

#include "function.h"

namespace housepi
{
	struct House
	{
		/// Household contact rate matrix.
		struct Contact
		{
			std::vector<std::vector<Attribute>> mtx = {}; /// contact rates
		};
	};
}
