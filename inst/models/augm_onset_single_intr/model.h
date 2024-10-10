/*
 * model.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////  * *** * AUGMENTED ONSET - SINGLE INTRODUCTION * *** *  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MODEL__H
#define MODEL__H

#include "profiles.h"
#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Model : class //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
public:
	//===== Model dependent methods =====//
	double transmission_rate(const Household &this_house) const;
	double foi(const Household &this_house, int i, const std::vector<squarematrix> &profiles) const;

	//===== Override General Methods =====//
	void init_model_for_inference() override;
	void init_model_for_simulation() override { init_model_for_inference(); };
	void update_augm_dependecies(Household &this_house, int i) override;

	//----- Rate Choices -----//
	void assign_group_rate_to_individual(Individual &this_indiv) override;
	void assign_contact_rate_to_household(Household &this_house, int from, int to) override;

	//----- Override Methods for Inference -----//
	double compute_house_loglike(const Household &this_house) const override;

	//----- Override Methods for Simulations -----//
	void simulate_single_house(Household &this_house) override;
};

#endif // MODEL__H