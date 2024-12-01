/*
 * model.h
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  * *** * E.G.AUGMENTED ONSET - SINGLE INTRODUCTION * *** *  ///////////////////
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
	double _houseROI_(const Household &house, int i, const squarematrix &profiles) const;
	double ROI_inst(const Household &house, int i) const;
	double ROI_cumul_infec(const Household &house, int i) const;
	double ROI_cumul_suscept(const Household &house, int i) const;
	double time_end_followup(const Household &house) const;

	//===== Override General Methods =====//
	void init_model_for_inference() override;
	void init_model_for_simulation() override;
	void update_augm_dependecies(Household &house, int i) override;
	double compute_house_loglike(Household &house) const override;
	void simulate_single_house(Household &house) override;

	//----- Rate Choices -----//
	void assign_group_rate_to_individual(Household &house, int i) override;
	void assign_contact_rate_to_household(Household &house, int from, int to) override;
};

#endif // HOUSEPI__MODEL__H