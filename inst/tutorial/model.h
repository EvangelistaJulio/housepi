#ifndef HOUSEPI__MODEL__H
#define HOUSEPI__MODEL__H

#include "mcmc.h"

class Model : public MCMC
{
public:
	// override general MCMC methods
	double compute_house_obs_log_lik(Household &house) override final;
	void simulate_single_house(Household &house) override final;
};

#endif