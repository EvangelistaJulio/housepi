// household.hpp //

#ifndef DEF_HOUSEHOLD__H
#define DEF_HOUSEHOLD__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "varia.hpp"

// model
#define HEALTHY 0
#define SYMPTOMATIC 1

/////////////////////////////////////////////////////////
///////// Individual : Class declaration ///////////////
////////////////////////////////////////////////////////

class Individual
{
private:
	std::unordered_map<std::string, Augmented> i_augm = {};
	int i_status = 0;

public:
	// Attributes (data)
	int i_indid = 0;
	int i_age = 0;
	int i_sex = 0;
	bool i_isolation = false;
	bool i_recovered = false;
	int i_vacc_status = 0;
	timetype i_study_period = INFINITY;
	timetype i_vacc_delay = INFINITY;
	timetype i_onset_time = INFINITY;

	// augmented variables behavior list (model)
	std::vector<std::string> i_augm_list = {};
	std::vector<std::string> i_augm_with_model = {};
	std::vector<std::string> i_augm_to_be_inferred = {};

	// --- Constructors --- //
	Individual(
		int indid,
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

	// --- Methods (default) --- //

	void add_augm(std::string augmID,
				  double new_val,
				  std::string model_name = UNDEF,
				  kwargs_ext model_args_ext = {},
				  std::string strat_name = UNDEF,
				  kwargs strat_args = {},
				  double value_0 = INFINITY,
				  kwargs value_range = {});

	// Getter
	bool get_status() const { return i_status; }
	double get_augm(std::string augmID) const { return i_augm.at(augmID).get_value(); }
	bool augm_infer(std::string augmID) { return i_augm.at(augmID).get_infer(); }

	// Setter
	void set_status(int status) { i_status = status; }
	void set_augm(std::string augmID, double new_val) { i_augm.at(augmID).set_value(new_val); }
	void init_augm(std::string augmID, std::mt19937_64 &rng) { set_augm(augmID, augm_model_rand(augmID, rng)); }

	// model distribution
	double augm_model_rand(std::string augmID, std::mt19937_64 &rng) { return i_augm.at(augmID).model_rand(rng); }
	double augm_model_log_pdf(std::string augmID);
	double augm_model_log_ratio(std::string augmID, double old_value, double new_val);

	// proposer
	double augm_propose(std::string augmID, std::mt19937_64 &rng);
	double augm_propose_log_ratio(std::string augmID, double old_value, double new_val);
	bool augm_valid(std::string augmID, double new_val) { return i_augm.at(augmID).valid(new_val); }
};

/////////////////////////////////////////////////////////
///////// Household : Class declaration ///////////////
////////////////////////////////////////////////////////

class Household
{
public:
	std::vector<Individual> h_residents;
	std::vector<std::vector<double>> h_cumul_profile = {}; // h_cumul_profile[i][j] : from i to j
	std::vector<std::vector<double>> h_inst_profile = {};

	// Attributes
	int house_index = -1;
	size_t h_size = 0;
	size_t h_first_case = -1;
	timetype h_infec_time = INFINITY;
	std::vector<size_t> h_healthy_list = {};
	std::vector<size_t> h_symptomatic_list = {};
	double h_loglike = 0;

	// Constructors
	Household() {}
	~Household(){};

	// --- Methods (default) --- //

	// Getter
	int get_index() const { return house_index; }
	int get_index(size_t idx) const { return h_residents[idx].i_indid; } // FIXME : to be gen with data
	int get_status(size_t idx) const { return h_residents[idx].get_status(); }
	std::vector<std::string> get_augm_list(size_t idx) const { return h_residents[idx].i_augm_list; }
	std::vector<std::string> get_augm_with_model(size_t idx) const { return h_residents[idx].i_augm_with_model; }
	std::vector<std::string> get_augm_to_be_inferred(size_t idx) const { return h_residents[idx].i_augm_to_be_inferred; }

	// Setter
	void set_index(int index) { house_index = index; }
	void set_status(size_t idx, int status);
	void add_individual(Individual this_indiv, int status);

	// Aliases (augmented variables)
	double get_augm(std::string augmID, size_t idx) const { return h_residents[idx].get_augm(augmID); }
	void set_augm(std::string augmID, size_t idx, double new_val) { h_residents[idx].set_augm(augmID, new_val); }

	// proposer
	double augm_propose(std::string augmID, size_t idx, std::mt19937_64 &rng) { return h_residents[idx].augm_propose(augmID, rng); }
	bool augm_valid(std::string augmID, size_t idx, double new_val) { return h_residents[idx].augm_valid(augmID, new_val); }
	double augm_propose_log_ratio(std::string augmID, size_t idx, double old_value, double new_val) { return h_residents[idx].augm_propose_log_ratio(augmID, old_value, new_val); }

	// model
	double augm_model_log_pdf(std::string augmID, size_t idx) { return h_residents[idx].augm_model_log_pdf(augmID); }
	double augm_model_log_ratio(std::string augmID, size_t idx, double old_value, double new_val) { return h_residents[idx].augm_model_log_ratio(augmID, old_value, new_val); }
};

#endif // DEF_HOUSEHOLD__H