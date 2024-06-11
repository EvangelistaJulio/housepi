// user_model.h //

#ifndef MODEL__H
#define MODEL__H

#include "mcmc.h"

// auxiliary types to optimize memory
typedef unsigned short usi;
typedef std::vector<usi> vec_usi;

//=====AUXILIARY FUNCTIONS DECLARATION=====//
// combinations look up table
unsigned binomial_coefficient(const usi N, const usi k);
void _recursive_combinations_span(usi N, usi k, vec_usi &current, usi start, usi loop, std::vector<vec_usi> &results);
void store_combinations(const usi N, const usi k);

// recursive likelihood formulation
double _recursive_function(usi N, std::vector<double> &scape_matrix);
double ll_saturation_heterogeneous(usi N, usi k, usi *indidx, std::vector<double> &escape_matrix);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
	void initialize_model() override;
	double compute_single_house_loglike(int h) override;
	void simulate_single_house(int h) override;
};

#endif // MODEL__H