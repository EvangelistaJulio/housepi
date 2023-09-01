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
extern Onset_to_Transm Infectivity_Profile;

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
	timetype i_study_period = INFINITY;
	timetype i_vacc_delay = INFINITY;
	timetype i_onset_time = INFINITY;
	Onset_to_Transm *ptr_to_infec_profile;

	std::vector<std::string> i_augm_list = {};
	std::vector<std::string> i_augm_with_model_list = {};

	// Constructors
	Individual(
		int indid,
		int status,
		int age,
		int sex,
		bool isolation,
		bool recovered,
		int vacc_status,
		timetype study_period,
		timetype vacc_delay,
		timetype onset_time);
	Individual() = default;
	~Individual(){};

	// --- Methods (default) ---

	void add_augm(std::string augmID,
				  double new_value,
				  std::string model_name = UNDEF,
				  kwargs_ext model_args_ext = {},
				  std::string strat_name = UNDEF,
				  kwargs strat_args = {},
				  kwargs value_range = {});

	// Getter
	double get_augm(std::string augmID) const { return i_augm.at(augmID).get_value(); }

	// Setter
	void set_augm(std::string augmID, double new_value) { i_augm.at(augmID).set_value(new_value); }

	// Computations
	double augm_model_rand(std::string augmID, std::mt19937_64 &rng) { return i_augm.at(augmID).model_rand(rng); }
	double augm_model_log_pdf(std::string augmID);
	double augm_model_log_ratio(std::string augmID, double new_value); // FIXME: old value always == value_0
	double augm_propose(std::string augmID, std::mt19937_64 &rng);
	double augm_propose_log_ratio(std::string augmID);
	bool augm_valid(std::string augmID, double new_value) { return i_augm.at(augmID)._base_.valid(new_value); }

	// --- Methods (model) ---
	timetype get_infec_time() const;
	void init_augm_data(std::mt19937_64 &rng);

	// Infectivity Profiles
	double infectivity_inst(timetype t_transm);
	double infectivity_cumul(timetype t_transm);
};

/////////////////////////////////////////////////////////
///////// Household : Class declaration ///////////////
////////////////////////////////////////////////////////

class Household
{
private:
	std::vector<Individual> h_residents;
	std::vector<std::vector<double>> h_cumul_profile = {}; // h_cumul_profile[i][j] : from i to j
	std::vector<std::vector<double>> h_inst_profile = {};

public:
	// Attributes
	size_t h_size = 0;
	size_t h_first_case = -1;
	timetype h_infec_time = INFINITY;
	std::vector<size_t> h_healthy_list = {};
	std::vector<size_t> h_symptomatic_list = {};
	Onset_to_Transm *ptr_to_infec_profile;
	double h_loglike = 0.0;

	// Constructors
	Household(){};
	~Household(){};

	// --- Methods (default) ---
	// Getter
	int get_status(size_t idx) const { return h_residents[idx].i_status; }
	std::vector<std::string> get_augm_list(size_t idx) const { return h_residents[idx].i_augm_list; }
	timetype get_infec_time(size_t idx) const { return h_residents[idx].get_infec_time(); }

	// Setter
	void add_individual(Individual this_indiv);
	void init_augm_data(std::mt19937_64 &rng);
	void init_household_status();
	void update_household_first_case();

	// Aliases (augmented variables)
	double get_augm(std::string augmID, size_t idx) const { return h_residents[idx].get_augm(augmID); }
	void set_augm(std::string augmID, size_t idx, double infec_time);
	double augm_propose(std::string augmID, size_t idx, std::mt19937_64 &rng) { return h_residents[idx].augm_propose(augmID, rng); }
	bool augm_valid(std::string augmID, size_t idx, double new_value) { return h_residents[idx].augm_valid(augmID, new_value); }

	double augm_model_log_pdf(std::string augmID, size_t idx) { return h_residents[idx].augm_model_log_pdf(augmID); }
	double augm_model_log_ratio(std::string augmID, size_t idx, double new_value) { return h_residents[idx].augm_model_log_ratio(augmID, new_value); }
	double augm_propose_log_ratio(std::string augmID, size_t idx) { return h_residents[idx].augm_propose_log_ratio(augmID); }

	// Computations (default)
	void init_infectivity_profiles(Onset_to_Transm *infectivity_profile);
	void update_infectivity_profiles(size_t ind);

	// --- Methods (model) ---
	double log_infec_proba(size_t idx);
	double log_surv_proba(size_t idx);
	double compute_loglike();

	double model_contact_modulation(size_t idx, size_t infectee);
	double model_susceptibility(size_t idx);
	double model_community_risk_inst(timetype t1);
	double model_community_risk_cumul(timetype t0, timetype t1);
	double model_transmissability();
	double model_risk_beta(const std::vector<std::vector<double>> &rate_vector, size_t idx);
};

#endif // DEF_HOUSEHOLD__H