// household.hpp //

#ifndef DEF_HOUSEHOLD__H
#define DEF_HOUSEHOLD__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "individual.hpp"

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
	std::vector<size_t> h_symptomatic_list = {};

	// Constructors
	Household(){};
	~Household(){};

	// Getter
	int get_status(size_t idx) const { return h_residents[idx].i_status; }
	double get_augm_incub_period(size_t idx) const { return h_residents[idx].get_augm_incub_period(); }
	double get_augm_infec_time(size_t idx) const { return h_residents[idx].get_augm_infec_time(); }
	double get_augm_onset_time(size_t idx) const { return h_residents[idx].get_augm_onset_time(); }

	// Setter
	void add_individual(Individual this_indiv);
	void init_augm_data(std::mt19937_64 &rand_gen);
	void update_household_status();

	// Computations (augmented)
	double propose_augm_incub_period(size_t idx, std::mt19937_64 &rand_gen);
	double propose_augm_onset_time(size_t idx, std::mt19937_64 &rand_gen);
	void set_augm_incub_period(size_t idx, double infec_time);
	void set_augm_onset_time(size_t idx, double onset_time);

	// Computations (default)
	void init_infectivity_profiles();
	void update_infectivity_profiles(size_t ind);

	double log_incub_period_pdf(size_t idx);
	double log_infec_proba(size_t idx);
	double log_surv_proba(size_t idx, double t0);
	double compute_loglike();

	// Model // NOTE: this should be external from the class Household
	double model_contact_modulation(size_t idx, size_t infectee);
	double model_susceptibility(size_t idx);
	double model_community_risk_inst(double t1);
	double model_community_risk_cumul(double t0, double t1);
	double model_transmissability();

};

#endif // DEF_HOUSEHOLD__H