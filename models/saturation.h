// saturation.h //

#ifndef MODEL__H
#define MODEL__H

#include "mcmc.h"

// auxiliary types
typedef std::vector<int> gentype;

// auxiliary functions
unsigned long long factorial(int n);
std::vector<gentype> saturation_count_span(int n_infec);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
	// Override MCMC methods
	void initialize_model() override;
	double compute_single_house_loglike(int h) override;
	void simulate_infection() override;
};

#endif // MODEL__H