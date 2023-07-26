#ifndef DEF_HOUSEHOLD__H
#define DEF_HOUSEHOLD__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "parameter.hpp"

#define VACC_ONSET 90.

//----------------------------------------------------------------------
// Individual : Class declaration
//----------------------------------------------------------------------

class Individual {
	public:
		// Attributes
		int i_indid = 0;
		int i_is_case = 0;
		int i_age = 0;
		int i_sex = 0;
		bool i_isolation = false;
		bool i_recovered = false;
		int i_vacc_status = 0;
		double i_start_follow_up = 0.0;
		double i_study_period = 0.0;
		double i_vacc_delay = 0.0;
		double i_onset_time = INFINITY;
		double i_augm_onset_time = INFINITY; 
		double i_augm_infec_time = INFINITY;

		// Constructors
		Individual(
			int indid,
			int is_case,
			int age,
			int sex,
			bool isolation,
			bool recovered,
			int vacc_status,
			double start_follow_up,
			double study_period,
			double vacc_delay,
			double onset_time
			);
		Individual() = default;
		~Individual() {};

		// Methods
		double infec_profile(double tinf, double study_period);
		double infec_cumulative(double tinf, double study_period);
};

//----------------------------------------------------------------------
// Household : Class declaration
//----------------------------------------------------------------------

class Household {
	// Attributes
	size_t h_size = 0;
	bool h_contact_pattern = false;
	double h_infec_time = INFINITY;

	std::vector<int> h_confCase;
	std::vector<Individual> h_residents;
	std::vector<double> h_inst_community_infec;
	std::vector<double> h_cumul_community_infec;
	std::vector<std::vector<double>> h_cumul_Lambda;
	std::vector<std::vector<double>> h_inst_Lambda;

	public:
		// Constructors
		Household(bool contact_pattern, std::vector<double> instval, std::vector<double>cumval);
		Household() = default;
		~Household() {};
		
		// Getter
		size_t get_size() const;
		size_t num_infec_cases() const;
		size_t num_symptomatic() const;
		std::string get_resident_index() const;
		std::vector<int> get_infec_list() const;
		double get_augm_onset_time(size_t idx) const;

		// Setter
		void add_individual(Individual this_indiv);
		void init_infec_time(size_t idx, std::mt19937_64& rand_gen, bool display = false);
		void init_onset_time(size_t idx, std::mt19937_64& rand_gen, bool display = false);

		// Computations
		void set_augm_infec_time(size_t idx, double infec_time, bool display = false);
		void set_augm_onset_time(size_t idx, double onset_time, bool display = false);
		void init_infec_profiles(bool display = false);
		void update_infec_profiles(size_t ind, bool display = false);
		double new_infec_time(size_t idx, std::mt19937_64& rand_gen);
		double new_onset_time(size_t idx, std::mt19937_64& rand_gen);

		double child_infectivity (size_t curr, size_t infectee);
		double model_susceptibility (size_t curr);
		double log_incub_proba(size_t idx, bool display = false);
		double compute_loglike(bool display = false);
		double log_infec_proba(size_t curr, bool display = false);
		double log_surv_proba(size_t curr, double t0, bool display = false);
};

#endif // DEF_HOUSEHOLD__H