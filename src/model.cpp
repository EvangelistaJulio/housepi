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
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * See the COPYRIGHT file for details.
 * 
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the COPYING file for details.
 *
 * model.h
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  * *** * BLANK MODEL * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/messages.h"
#include "model.h"

/// Raise error when "blank" model.
void Model::init_model_for_inference() { hsp::msg::warning(
	"Computed log-likelihood with 'blank' model."); }

/// Raise error when "blank" model.
void Model::init_model_for_simulation() { hsp::msg::warning(
	"Performed simulation with 'blank' model."); }

/// Unitary likelihood.
double Model::compute_household_obs_log_lik(Household &) { return 0.; }

/// Augmented variable initialization.
void Model::simulate_single_household(Household &house)
{
	for (Individual &member : house.infected())
		init_missing_augm(member);
}
