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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Model : class //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
public:
	//===== Global Variables =====//
	//----- Infectivity Profile Functions -----//
	Infec_Prof_Skew_Logistic *infec_prof_ptr;

	//----- Attributes for Look-Up-Tables -----//
	vector<squarematrix> LUT_inst_profile = {};
	vector<squarematrix> LUT_cumul_profile = {};
	vector<int> LUT_house_first_case = {};

	//===== Model dependent methods =====//
	double transmission_rate(int h) const;
	double force_of_infec(int h, int i, const vector<squarematrix> &profiles) const;

	//===== Override General Methods =====//
	void init_model_for_inference() override;
	void init_model_for_simulation() override { init_model_for_inference(); };
	
	void augm_store_dependecies(int h, int i) override;

	//----- Rate Choices -----//
	void assign_rate_per_group_to_individual(Individual &this_indiv) override;
	void assign_rate_per_contact_to_household(int h, int from, int to) override;

	//----- Override Methods for Inference -----//
	double compute_house_loglike(int h) const override;

	//----- Override Methods for Simulations -----//
	void simulate_single_house(int h) override;
};

#endif // MODEL__H