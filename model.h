// model.h //

#ifndef MODEL__H
#define MODEL__H

#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "individual.h"

extern bool Verbose;

////////////////////////////////////////////////////////////////////////
// ////////////////// Extra files for the model //////////////////// //
///////////////////////////////////////////////////////////////////////

// Shared Parameter Dictionary // UPDATE : generalize for external files
void build_extra(std::string data_pathtofile);

// forward declaration of the Infectivity Profile Function
class Infectivity_Profile_Function;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public Parameter_Dictionary
{
protected:
	//=====Required Attributes=====//
	Individual *indiv;
	std::size_t num_indiv = 0;

	Infectivity_Profile_Function *infec_profile_func;
	std::set<std::string> augm_set = {};
	std::set<std::string> augm_with_model_set = {};
	std::set<std::string> augm_to_infer_set = {};
	std::vector<std::string> group_rate_list = {};
	std::vector<std::string> data_headers = {};
	
	//=====Optional Attributes=====//
	Household *househ;
	std::size_t num_househ = 0;

	double *inst_profile_lookuptab;
	double *cumul_profile_lookuptab;

	const double **contact_patterns;

public:
	//=====Required Methods=====//
	int get_status(std::size_t i) const { return indiv[i].get_status(); }
	
	void assign_infec_profile_function();
	//virtual void assign_augmented_variables(std::mt19937_64 &rng); // UPDATE: temporary
	void assign_rate_per_group_to_individual(Individual *this_indiv);
	void fill_model_lists(); // UPDATE: temporary

	void store_augm_var_dependecies(std::size_t i);

	double compute_loglike(std::size_t h);

	void simulate_infection(double delta_t, std::mt19937_64 &rng);

	//=====Optional Methods=====//

	std::size_t from_to(std::size_t i, std::size_t j) { return i * num_indiv + j; }
	void assign_rate_per_contact(std::size_t i, std::size_t k);
	double get_contact_from_to(std::size_t i, std::size_t k) { return (*contact_patterns[from_to(i, k)]); }

	//=====Infectivity Profiles=====//
	double compute_infec_profile_inst(std::size_t i, double t_transm);
	double compute_infec_profile_cumul(std::size_t i, double t_transm);

	//=====LogLikelihood=====//
	double _sum_hazards_(double *profiles, std::size_t idx);
};

#endif // MODEL__H