/*
 * model.h
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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
	double _houseROI_(const Household &house, indiv i, const matrix<double> &profiles) const;
	double ROI_inst(const Household &house, indiv i) const;
	double ROI_cumul_infec(const Household &house, indiv i) const;
	double ROI_cumul_suscept(const Household &house, indiv i) const;
	double time_end_followup(const Household &house) const;

	//===== Override General Methods =====//
	void init_model_for_inference() override;
	void init_model_for_simulation() override;
	void update_augm_dependecies(Household &house, indiv i) override;
	double compute_house_loglike(Household &house) const override;
	void simulate_single_house(Household &house) override;

	//----- Rate Choices -----//
	void assign_group_rate_to_individual(Household &house, indiv i) override;
	void assign_contact_rate_to_household(Household &house, indiv from, indiv to) override;
};

#endif // HOUSEPI__MODEL__H