#ifndef DEF_HOUSEHOLD__H
#define DEF_HOUSEHOLD__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "parameter.hpp"

#define NOT_A_CASE 0
#define SYMPTOMATIC 1
#define NOT_SYMPTOM 2
#define SYMPTOM_UNK_ONSET 3
#define VACC_ONSET 90.

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Community_Infec;
extern std::vector<double> Inst_Community_Infec;

//----------------------------------------------------------------------
// Individual : Class declaration
//----------------------------------------------------------------------

class Individual
{
	double i_augm_onset_time = INFINITY;
	double i_augm_infec_time = INFINITY;

public:
	// Attributes
	int i_indid = 0;
	int i_status = 0;
	int i_age = 0;
	int i_sex = 0;
	bool i_isolation = false;
	bool i_recovered = false;
	int i_vacc_status = 0;
	double i_study_period = 0.0;
	double i_vacc_delay = 0.0;
	double i_onset_time = INFINITY;

	// Constructors
	Individual(
		int indid,
		int status,
		int age,
		int sex,
		bool isolation,
		bool recovered,
		int vacc_status,
		double study_period,
		double vacc_delay,
		double onset_time);
	Individual() = default;
	~Individual(){};

	// Getter
	double get_augm_infec_time() const { return i_augm_infec_time; }
	double get_augm_onset_time() const { return i_augm_onset_time; }

	// Setter
	void init_augm_data(std::mt19937_64 &rand_gen);
	void set_augm_infec_time(double new_infec_time) { i_augm_infec_time = new_infec_time; }
	void set_augm_onset_time(double new_onset_time) { i_augm_onset_time = new_onset_time; }

	// Computations
	double propose_augm_infec_time(std::mt19937_64 &rand_gen);
	double propose_augm_onset_time(std::mt19937_64 &rand_gen);
	double infec_inst(double tinf, double study_period);
	double infec_cumul(double tinf, double study_period);
};

//----------------------------------------------------------------------
// Household : Class declaration
//----------------------------------------------------------------------

class Household
{
	std::vector<Individual> h_residents;
	std::vector<std::vector<double>> h_cumul_rate = {}; // h_cumul_rate[i][j] : from i to j
	std::vector<std::vector<double>> h_inst_rate = {};

public:
	// Attributes
	size_t h_size = 0;
	size_t h_first_case = -1;
	double h_infec_time = INFINITY;

	// Constructors
	Household(){};
	~Household(){};

	// Getter
	double get_augm_infec_time(size_t idx) const { return h_residents[idx].get_augm_infec_time(); }
	double get_augm_onset_time(size_t idx) const { return h_residents[idx].get_augm_onset_time(); }

	// Setter
	void add_individual(Individual this_indiv);
	void model_init_augm_data(std::mt19937_64 &rand_gen);

	// Computations
	void set_augm_infec_time(size_t idx, double infec_time);
	void set_augm_onset_time(size_t idx, double onset_time);
	void init_infec_profiles();
	void update_h_infec_time();
	void update_infec_profiles(size_t ind);
	double propose_augm_infec_time(size_t idx, std::mt19937_64 &rand_gen);
	double propose_augm_onset_time(size_t idx, std::mt19937_64 &rand_gen);

	// Model // NOTE: this should be external from the class Household
	double model_contact_modulation(size_t idx, size_t infectee);
	double model_susceptibility(size_t idx);
	double model_community_risk_inst(double t1);
	double model_community_risk_cumul(double t0, double t1);
	double model_transmissability();
	double log_incub_proba(size_t idx);
	double compute_loglike();
	double log_infec_proba(size_t idx);
	double log_surv_proba(size_t idx, double t0);
};

#endif // DEF_HOUSEHOLD__H