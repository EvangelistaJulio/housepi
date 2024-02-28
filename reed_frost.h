// reed_frost.h //

#ifndef MODEL__H
#define MODEL__H

#include <string>
#include <cmath>
#include <random>
#include <vector>
#include <map>
#include <iostream>

#include "mcmc.h"

typedef std::vector<int> gentype;

unsigned long long factorial(int n);
std::vector<gentype> generation_span(int n_infec);
std::vector<gentype> generation_count_span(int n_infec);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
protected:
	std::map<int, std::vector<gentype>> generations_Look_Up_Table = {};

public:
	// Constructors
	Model(Setup_MCMC_Struct new_setup) : MCMC(new_setup) {
	for (int k = 1; k < 10; k++)
		{
			generations_Look_Up_Table.emplace(k, generation_count_span(k));
		}
	};
	~Model() {}

	double compute_single_house_loglike(std::size_t h) override;
	void simulate_infection() override;
};

#endif // MODEL__H