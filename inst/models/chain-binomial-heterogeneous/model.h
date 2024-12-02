/*
 * model.h
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef HOUSEPI__MODEL__H
#define HOUSEPI__MODEL__H

#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Model : class //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
	//===== Model dependent methods =====//
	double community_risk() const;
	double person_to_person_risk() const;
	void fill_house_escape_matrix(Household &house) const;

	//===== Override MCMC methods =====//
	void init_model_for_inference() override;
	void init_model_for_simulation() override;
	double compute_house_loglike(Household &house) const override;
	void simulate_single_house(Household &house) override;
	void assign_group_rate_to_individual(Household &house, indiv i) override;
};

#endif // HOUSEPI__MODEL__H