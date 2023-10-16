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
	Onset_to_Transm(kwargs new_args) : Base_Distrib(new_args) {}
	Onset_to_Transm(kwargs_ext new_args_ext) : Base_Distrib(new_args_ext) {}
	~Onset_to_Transm(){};

	double pdf(timetype tost, timetype t_incub);
	double cdf(timetype tost, timetype t_incub);
};

///////////////////////////////////////////////////////////////////////
// //////////////// Modified Individual : class //////////////////// //
///////////////////////////////////////////////////////////////////////

class Individual : public Base_Individual
{
public:
	size_t home = -1;
	bool is_child = false;
	bool is_adult = false;
	bool is_isol = false;
	bool is_vacc_0 = false;
	bool is_vacc_1 = false;
	bool is_vacc_2 = false;

	using Base_Individual::Base_Individual;
};

///////////////////////////////////////////////////////////////////////
// ///////////////////////// Model : class ///////////////////////// //
///////////////////////////////////////////////////////////////////////

class Model : public Parameter_Dictionary
{
protected:
	Individual *indiv;
	Household *househ;

	// infectivity profile
	Onset_to_Transm *infec_profile_func;
	double *inst_profile;
	double *cumul_profile;

	// augmented variables behavior lists
	std::vector<std::string> augm_list = {};
	std::vector<std::string> augm_with_model = {};
	std::vector<std::string> augm_to_be_inferred = {};

public:
	size_t num_indiv = 0;
	size_t num_symptomatic = 0;
	size_t num_househ = 0;

	std::vector<std::string> data_headers = {};

	//=====Methods=====//
	void pre_process_data();
	int get_status(size_t i) const { return indiv[i].get_status(); }

	timetype study_time_start();
	timetype study_time_end();

	//=====Augmetend Data=====//
	void assign_augm_data(std::mt19937_64 &rng);
	void _assign_augm_to_ind_(size_t idx, std::mt19937_64 &rng); // UPDATE: temporary
	timetype get_infec_time(size_t i) const;

	//=====Infectivity Profile=====//
	size_t from_to(int i, int j) { return i * num_indiv + j; }
	void initialize_infectivity_profile();
	double infectivity_inst(size_t i, timetype t_transm);
	double infectivity_cumul(size_t i, timetype t_transm);
	void update_infectivity_profiles(size_t i);

	void update_first_case(size_t h);

	//=====LogLikelihood=====//
	double risk_beta(double *rate_vector, size_t idx);
	double risk_surv(size_t idx);
	double risk_infec(size_t idx);
	double compute_loglike(size_t h);

	double community_risk_inst(timetype t1);
	double community_risk_cumul(timetype t0, timetype t1);
	double contact_modulation(size_t infectee, size_t idx);
	double susceptibility(size_t idx);
	double transmissability(size_t i);

	void simulate_infection(timetype delta_t, std::mt19937_64 &rng);
};

#endif // MODEL__H