/*
 * model.h
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////  * *** * AUGMENT INFECTION TIMES - SINGLE INTRODUCTION * *** *  /////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HOUSEPI__MODEL__H
#define HOUSEPI__MODEL__H

#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Model : class //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
public:
	//===== Model dependent methods =====//
	double transmission_rate(const Household &house) const;

	double loglike_healthy(const Household &house) const;
	double loglike_infectious(const Household &house) const;
	double roi_inst_simulation(const std::vector<size_t> &infectors,
							   Household &house,
							   Individual &member_i,
							   double t_infec_i) const;

	//===== Override General Methods =====//
	void init_user_advanced_setup() override final;
	void init_model_for_inference() override final;
	void init_model_for_simulation() override final;
	void update_cache_dependecies(Household &house, Individual &member_i) override final;
	double compute_household_obs_loglike(Household &house) override final;
	void simulate_single_household(Household &house) override final;
};

#endif // HOUSEPI__MODEL__H