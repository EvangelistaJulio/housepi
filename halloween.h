// halloween.h //

#ifndef HALLOW__H
#define MODEL__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "mcmc.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
protected:
	double *inst_profile_lookuptab;
	double *cumul_profile_lookuptab;

public:
	// Constructors
	Model(Setup_MCMC_Struct new_setup) : MCMC(new_setup) {};
	~Model() {}

	//=====Required Methods=====//
	
	void assign_infec_profile_function() override;
	void assign_rate_per_group_to_individual(Individual *this_indiv) override;
	void store_augm_var_dependecies(std::size_t i) override;
	void initialize_model() override;

	double compute_single_house_loglike(std::size_t h) override;

	void simulate_infection() override;

	//=====Optional Methods=====//
	void assign_rate_per_contact(std::size_t i, std::size_t k) override;

	//=====Infectivity Profiles=====//
	double compute_infec_profile_inst(std::size_t i, double t_transm);
	double compute_infec_profile_cumul(std::size_t i, double t_transm);

	//=====LogLikelihood=====//
	double _sum_hazards_(double *profiles, std::size_t idx);
};

#endif // MODEL__H