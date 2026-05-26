/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>

#include "tools/status.h"
#include "tools/loglike.h"

#include "function.h"
#include "variable.h"

namespace housepi
{
	struct Indiv
	{
		/// Individual wide-format data.
		struct Wide
		{
			int status = status::MISSING;
			std::vector<double> vec = {};
		};

		/// Individual long-format data.
		struct Long
		{
			std::vector<std::vector<double>> vec = {};
		};

		/// Individual augmented variables.
		struct Augm
		{
			size_t idx = 0;					/// index of status specific augmented definition
			std::vector<Variable> vec = {}; /// augmented variable definitions
			LogLike log_lik;				/// augmented model log-likelihood
		};

		/// Individual group rates.
		struct Group
		{
			std::vector<Attribute> vec = {}; /// group rates
		};
	};
}
