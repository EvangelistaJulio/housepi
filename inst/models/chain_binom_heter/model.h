// model.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MODEL__H
#define MODEL__H

#include "binomial.h"
#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////// Model : class //////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
	//===== Model dependent methods =====//
	void fill_house_escape_matrix(const Household &this_house, std::vector<double> &commu_esc_proba, squarematrix &p2p_esc_proba) const;

	//===== Override MCMC methods =====//
	void init_model_for_inference() override;
	double compute_house_loglike(const Household &this_house) const override;
	void simulate_single_house(Household &this_house) override;
	void assign_group_rate_to_individual(Individual &this_indiv) override;
	void assign_contact_rate_to_household(Household &this_house, int i, int k) override;
};

#endif // MODEL__H