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
private:
	//----- Infectivity Profile Functions -----//
	std::shared_ptr<Infec_Prof_Ferretti20> infec_prof_ptr;

	//----- Attributes for Look-Up-Tables -----//
	std::vector<squarematrix> LUT_inst_prof;
	std::vector<squarematrix> LUT_cumul_prof;
	std::vector<int> LUT_hh_first_case;

public:
	//===== Model dependent methods =====//
	double transmission_rate(const Household &house) const;
	double foi(const Household &house, int i, const std::vector<squarematrix> &profiles) const;

	//===== Override General Methods =====//
	void init_model_for_inference() override;
	void init_model_for_simulation() override { init_model_for_inference(); };
	void update_augm_dependecies(Household &house, int i) override;

	//----- Rate Choices -----//
	void assign_group_rate_to_individual(Household &house, int i) override;
	void assign_contact_rate_to_household(Household &house, int from, int to) override;

	//----- Override Methods for Inference -----//
	double compute_house_loglike(const Household &house) const override;

	//----- Override Methods for Simulations -----//
	void simulate_single_house(Household &house) override;
};

#endif // MODEL__H