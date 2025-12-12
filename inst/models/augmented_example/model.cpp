/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Author: Francesco Camaglia
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 * See the COPYRIGHT file for details.
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the COPYING file for details.
 *
 * model.cpp
 * Model template of the 'housepi' package.
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  * *** * DATA AUGMENTATION EXAMPLE * *** *  ///////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"
// FLAG
#include "tools/messages.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute the log-likelihood within household.
double Model::compute_household_obs_log_lik(Household &house)
{
	double output = 0.0;
	for (const Individual &member : house.infected())
	{
		// FLAG
		hsp::msg::info(member.id() + ": " +
					   member.get_augm_variable("proba").str_summary() + "\n");
		output += std::log(member.augm("proba"));
	}

	return output;
}

/// Augmented variable initialization.
void Model::simulate_single_household(Household &house)
{
	for (Individual &member : house.infected())
		init_missing_augm(member);
}
