/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <vector>

#include "tools/status.h"
#include "tools/loglike.h"

#include "function.h"
#include "variable.h"

namespace hsp
{
	struct Indiv
	{
		/// Individual wide-format data.
		struct Wide
		{
			int status = status::MISSING;
			std::vector<double> data_vec = {};
		};

		/// Individual long-format data.
		struct Long
		{
			std::vector<std::vector<double>> data_vec = {};
		};

		/// Individual augmented variables.
		struct Augm
		{
			size_t augm_idx = 0;				 /// index of status specific augmented definition
			std::vector<Variable> augm_vec = {}; /// augmented variable definitions
			LogLike augm_log_lik;				 /// augmented model log-likelihood
		};

		/// Individual group rates.
		struct Group
		{
			std::vector<Attribute> group_rate_vec = {}; /// group rates
		};
	};
}
