// model.hpp //

#ifndef MODEL__H
#define MODEL__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "individual.h"

///////////////////////////////////////////////////////////////////////
// //////////////// Onset to Transmission : class ////////////////// //
///////////////////////////////////////////////////////////////////////

class Onset_to_Transm : public Base_Distrib
{
public:
	// Constructors
	Onset_to_Transm(kwargs new_args) : Base_Distrib(new_args, {"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}) {}
	Onset_to_Transm(kwargs_ext new_args_ext) : Base_Distrib(new_args_ext, {"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}) {}
	~Onset_to_Transm(){};

	double pdf(timetype tost, timetype t_incub);
	double cdf(timetype tost, timetype t_incub);
};

///////////////////////////////////////////////////////////////////////
// ///////////////////////// Model : class ///////////////////////// //
///////////////////////////////////////////////////////////////////////

class Model : public Parameter_Dictionary
{
protected:
	Individual *indiv;
	Household *househ;

	// infectivity profiles
	Onset_to_Transm *infec_profile_func;
	double *inst_profile;
	double *cumul_profile;

	const double **contact_patterns;

	// augmented variables behavior lists
	std::vector<std::string> augm_list = {};
	std::vector<std::string> augm_with_model = {};
	std::vector<std::string> augm_to_be_inferred = {};

	std::vector<std::string> rate_list = {};

public:
	size_t num_indiv = 0;
	size_t num_househ = 0;

	std::vector<std::string> data_headers = {};

	//=====Methods=====//
	int get_status(size_t i) const { return indiv[i].get_status(); }
	double get_contact_from_to(size_t i, size_t k) { return (*contact_patterns[from_to(i, k)]); }

	//=====Augmetend Data=====//
	void assign_augm_to_individual_(Individual *this_indiv, std::mt19937_64 &rng); // UPDATE: temporary

	void assign_rate_to_individual(Individual *this_indiv);
	void assign_contact_patterns(size_t i, size_t k);

	//=====Structure Factors=====//
	double transmissability(size_t h);

	//=====Infectivity Profiles=====//
	size_t from_to(size_t i, size_t j) { return i * num_indiv + j; }
	void initialize_infec_profiles();
	double infec_profile_inst(size_t i, timetype t_transm);
	double infec_profile_cumul(size_t i, timetype t_transm);
	void update_infec_profile(size_t i);
	void fill_model_lists();

	//=====LogLikelihood=====//
	double _sum_hazards_(double *profiles, size_t idx);
	double risk_surv(size_t idx);
	double risk_infec(size_t idx);
	double compute_loglike(size_t h);

	double community_risk_inst(timetype t1);
	double community_risk_cumul(timetype t0, timetype t1);

	//=====Simulation=====//
	void simulate_infection(timetype delta_t, std::mt19937_64 &rng);
};

#endif // MODEL__H