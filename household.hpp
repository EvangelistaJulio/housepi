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

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Community_Infec;
extern std::vector<double> Inst_Community_Infec;

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

	std::vector<Individual> h_residents;
	std::vector<std::vector<double>> h_cumul_Lambda = {};
	std::vector<std::vector<double>> h_inst_Lambda = {};

	public:
		bool h_verbose = false;

		// Constructors
		Household(bool contact_pattern,	bool verbose);
		Household() = default;
		~Household() {};
		
		// Getter
		size_t get_size () const;
		size_t num_symptomatic () const;
		std::string get_resident_index () const;
		double get_augm_onset_time (size_t idx) const;

		// Setter
		void add_individual (Individual this_indiv);
		void init_infec_time (size_t idx, std::mt19937_64& rand_gen);
		void init_onset_time (size_t idx, std::mt19937_64& rand_gen);

		// Computations
		void set_augm_infec_time (size_t idx, double infec_time);
		void set_augm_onset_time (size_t idx, double onset_time);
		void init_infec_profiles ();
		void update_infec_profiles (size_t ind);
		double new_infec_time (size_t idx, std::mt19937_64& rand_gen);
		double new_onset_time (size_t idx, std::mt19937_64& rand_gen);

		double child_infectivity (size_t curr, size_t infectee);
		double model_susceptibility (size_t curr);
		double log_incub_proba (size_t idx);
		double compute_loglike ();
		double log_infec_proba (size_t curr);
		double log_surv_proba (size_t curr, double t0);
};

#endif // DEF_HOUSEHOLD__H