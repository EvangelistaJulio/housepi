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

#include <string>
#include <vector>
#include <iostream>

#include "profiles.h"
#include "mcmc.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

// alias for a square matrix
typedef std::vector<std::vector<double>> squarematrix;

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
	//===== Override General Methods =====//
	void common_init_model();
	void init_model_for_inference() override { common_init_model(); };
	void init_model_for_simulation() override { common_init_model(); };
	void augm_store_dependecies(int h, int i) override;
	double get_infec_time(int h, int i) const override;

	//----- Override Methods for Inference -----//
	double compute_house_loglike(int h) override;

	//----- Override Methods for Simulations -----//
	void simulate_single_house(int h) override;

	//----- Rate Choices -----//
	void assign_rate_per_group_to_individual(Individual &this_indiv) override;
	void assign_rate_per_contact_to_household(int h, int from, int to) override;
	
	//===== Model dependent methods =====//
	double compute_infec_prof_inst(int h, int i, double t_transm);
	double compute_infec_prof_cumul(int h, int i, double t_transm);
	double hazards_of_infec(int h, int i, std::vector<squarematrix> &profiles);
};

#endif // MODEL__H