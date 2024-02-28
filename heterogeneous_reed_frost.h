// heterogeneous_reed_frost.h //

#ifndef MODEL__H
#define MODEL__H

#include <iostream>
#include <string>
#include <cmath>
#include <random>
#include <vector>
#include <map>

#include "mcmc.h"

// auxiliary types to optimize memory
typedef unsigned short usi;
typedef std::vector<usi> vec_usi;
typedef std::vector<std::vector<double>> mtx;

// AUXILIARY FUNCTIONS DECLARATION
double Recursive_Function(usi N, mtx &escape_matrix);
double Likelihood_Heterogeneous_ReedFrost(usi N, usi k, vec_usi &indiv, mtx &escape_matrix);
unsigned binomial_coefficient(const usi N, const usi k);
void _recursive_combinations_span(usi N, usi k, vec_usi &current, usi start, usi loop, std::vector<vec_usi> &results);
void store_combinations(const usi N, const usi k);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
public:
	// Constructors
	Model(Setup_MCMC_Struct new_setup) : MCMC(new_setup) {};
	~Model() {}

	void initialize_model() override;
	double compute_single_house_loglike(std::size_t h) override;
	void simulate_infection() override;

};

#endif // MODEL__H