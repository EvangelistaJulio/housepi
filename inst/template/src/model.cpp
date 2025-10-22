/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * model.h
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  * *** * BLANK MODEL * *** *  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "tools/error.h"
#include "model.h"

/// Raise error when "blank" model.
void Model::init_model_for_inference() { hsp::msg::warning(
	"Computed loglikelihood with 'blank' model."); }

/// Raise error when "blank" model.
void Model::init_model_for_simulation() { hsp::msg::warning(
	"Performed simulation with 'blank' model."); }

/// Unitary likelihood.
double Model::compute_household_obs_loglike(Household &) { return 0.; }

/// Augmented variable initialization.
void Model::simulate_single_household(Household &house)
{
	for (Individual &member : house.infected())
		init_missing_augm(member);
}
