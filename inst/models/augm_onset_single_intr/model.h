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
	Infec_Prof_Ferretti20 *infec_prof_ptr;

	//----- Attributes for Look-Up-Tables -----//
	vector<squarematrix> LUT_inst_profile = {};
	vector<squarematrix> LUT_cumul_profile = {};
	vector<int> LUT_house_first_case = {};

	//===== Model dependent methods =====//
	double transmission_rate(const Household &this_house) const;
	double force_of_infec(const Household &this_house, int i, const vector<squarematrix> &profiles) const;

	double inst_FOI(const Household &this_house, int i) const { return force_of_infec(this_house, i, LUT_inst_profile); }
	double cumul_FOI(const Household &this_house, int i) const { return force_of_infec(this_house, i, LUT_cumul_profile); }

	//===== Override General Methods =====//
	void init_model_for_inference() override;
	void init_model_for_simulation() override { init_model_for_inference(); };
	void augm_store_dependecies(Household &this_house, int i) override;

	//----- Rate Choices -----//
	void assign_rate_per_group_to_individual(Individual &this_indiv) override;
	void assign_rate_per_contact_to_household(Household &this_house, int from, int to) override;

	//----- Override Methods for Inference -----//
	double compute_house_loglike(const Household &this_house) const override;

	//----- Override Methods for Simulations -----//
	void simulate_single_house(Household &this_house) override;
};

#endif // MODEL__H