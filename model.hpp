// model.hpp //

#ifndef MODEL__H
#define MODEL__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "individual.hpp"

///////////////////////////////////////////////////////////////////////
// //////////////// Modified Individual : class //////////////////// //
///////////////////////////////////////////////////////////////////////

class Individual : public Base_Individual
{
public:
	size_t home;
	bool is_child = false;
	bool is_mother = false;
	bool is_father = false;
	bool is_recov = false;
	bool is_isol = false;
	bool is_vacc_0 = false;
	bool is_vacc_1 = false;
	bool is_vacc_1_V_onset = false;
	bool is_vacc_2 = false;

	using Base_Individual::Base_Individual;
};

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
// ///////////////////////// Model : class ///////////////////////// //
///////////////////////////////////////////////////////////////////////

class Model : public Parameter_Dictionary
{
protected:
	Onset_to_Transm *infectivity_profile;
	std::vector<Individual> indiv;
	std::vector<Household> househ;

	// augmented variables behavior lists
	std::vector<std::string> augm_list = {};
	std::vector<std::string> augm_with_model = {};
	std::vector<std::string> augm_to_be_inferred = {};

	// infectivity profile : profile[i][j] : from i to j
	std::vector<std::vector<double>> cumul_profile = {}; 
	std::vector<std::vector<double>> inst_profile = {};

public:
	size_t num_indiv = 0;
	size_t num_symptomatic = 0;
	size_t num_househ = 0;

	void add_individual(Individual this_indiv);
	void assign_classes();
	void assign_community_structure();
	int get_status(size_t i) const { return indiv[i].get_status(); }

	timetype max_study_period();

	// Augmetend Data
	void assign_augm_data(std::mt19937_64 &rng);
	void _assign_augm_to_ind_(size_t idx, std::mt19937_64 &rng); // UPDATE: temporary
	timetype get_infec_time(size_t i) const;

	// Infectivity Profile
	void assign_infec_profile();
	double infectivity_inst(size_t i, timetype t_transm);
	double infectivity_cumul(size_t i, timetype t_transm);
	void clear_infectivity_profiles();
	void update_infectivity_profiles(size_t i);
	void update_first_case(size_t h);

	double risk_beta(std::vector<std::vector<double>> &rate_vector, size_t idx);
	double log_surv_proba(size_t idx);
	double log_infec_proba(size_t idx);
	double compute_loglike(size_t h);

	double community_risk_inst(timetype t1);
	double community_risk_cumul(timetype t0, timetype t1);
	double contact_modulation(size_t infectee, size_t idx);
	double susceptibility(size_t idx);
	double transmissability(size_t i);

	//void simulate_data(timetype delta_t, std::mt19937_64 &rng);
};

#endif // MODEL__H