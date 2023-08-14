// household.hpp //

#ifndef DEF_HOUSEHOLD__H
#define DEF_HOUSEHOLD__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "parameter.hpp"

// List of augmented variables
extern std::vector<std::string> List_Augm_Var;

/////////////////////////////////////////////////////////
///////// Individual : Class declaration ///////////////
////////////////////////////////////////////////////////

class Individual
{
private:
	std::unordered_map<std::string, Augmented> i_augm = {};

public:
	// Attributes
	int i_indid = 0;
	int i_status = 0;
	int i_age = 0;
	int i_sex = 0;
	bool i_isolation = false;
	bool i_recovered = false;
	int i_vacc_status = 0;
	double i_study_period = INFINITY;
	double i_vacc_delay = INFINITY;
	double i_onset_time = INFINITY;

	std::vector<std::string> i_augm_list = {};

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

	// --- Methods (default) ---
	// Getter
	double get_augm(std::string augmID) const { return i_augm.at(augmID).get_value(); }

	// Setter
	void set_augm(std::string augmID, double new_value) { i_augm.at(augmID).set_value(new_value); }

	// Computations
	bool valid_augm(std::string augmID, double new_value) { return i_augm.at(augmID).valid(new_value); }
	double propose_augm(std::string augmID, std::mt19937_64 &rand_gen);
	double log_ratio_prior_augm(std::string augmID, double new_value);
	double log_ratio_proposal_augm(std::string augmID, double new_value);

	// --- Methods (model) ---
	double get_infec_time() const;
	void init_augm_data(std::mt19937_64 &rand_gen);

	// Infectivity Profiles
	double infectivity_inst(double t_transm);
	double infectivity_cumul(double t_transm);
};

/////////////////////////////////////////////////////////
///////// Household : Class declaration ///////////////
////////////////////////////////////////////////////////

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
	std::vector<size_t> h_healthy_list = {};
	std::vector<size_t> h_symptomatic_list = {};

	// Constructors
	Household(){};
	~Household(){};

	// --- Methods (default) ---
	// Getter
	int get_status(size_t idx) const { return h_residents[idx].i_status; }
	std::vector<std::string> get_augm_list(size_t idx) const { return h_residents[idx].i_augm_list; }
	double get_infec_time(size_t idx) const { return h_residents[idx].get_infec_time(); }

	// Setter
	void add_individual(Individual this_indiv);
	void init_augm_data(std::mt19937_64 &rand_gen);
	void update_household_status();

	// Aliases (augmented variables)
	double get_augm(std::string augmID, size_t idx) const { return h_residents[idx].get_augm(augmID); }
	double propose_augm(std::string augmID, size_t idx, std::mt19937_64 &rand_gen);
	void set_augm(std::string augmID, size_t idx, double infec_time);
	bool valid_augm(std::string augmID, size_t idx, double new_value);
	double log_ratio_prior_augm(std::string augmID, size_t idx, double new_value) { return h_residents[idx].log_ratio_prior_augm(augmID, new_value); }
	double log_ratio_proposal_augm(std::string augmID, size_t idx, double new_value) { return h_residents[idx].log_ratio_proposal_augm(augmID, new_value); }

	// Computations (default)
	void init_infectivity_profiles();
	void update_infectivity_profiles(size_t ind);

	// --- Methods (model) ---
	double log_incub_proba(size_t idx); // FIXME : why not a augm_var prior?
	double log_infec_proba(size_t idx);
	double log_surv_proba(size_t idx);
	double compute_loglike();

	double model_contact_modulation(size_t idx, size_t infectee);
	double model_susceptibility(size_t idx);
	double model_community_risk_inst(double t1);
	double model_community_risk_cumul(double t0, double t1);
	double model_transmissability();
	double model_risk_beta(const std::vector<std::vector<double>> &rate_vector, size_t idx);
};

#endif // DEF_HOUSEHOLD__H