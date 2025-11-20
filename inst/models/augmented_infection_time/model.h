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
 * model.h
 * Model template of the 'housepi' package.
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////  * *** * AUGMENT INFECTION TIMES - SINGLE INTRODUCTION * *** *  /////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HOUSEPI__MODEL__H
#define HOUSEPI__MODEL__H

#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Model : class //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
public:
	// model dependent methods
	double transmission_rate(const Household &house) const;

	double loglike_healthy(const Household &house) const;
	double loglike_infectious(const Household &house) const;
	double roi_inst_simulation(Household &house,
							   Individual &member_i,
							   double t_infec_i) const;

	// override general MCMC methods
	void init_advanced_setup() override final;
	void init_model_for_inference() override final;
	void init_model_for_simulation() override final;
	void update_cache_dependecies(Household &house, Individual &member_i) override final;
	double compute_household_obs_log_lik(Household &house) override final;
	void simulate_single_household(Household &house) override final;
};

#endif // HOUSEPI__MODEL__H