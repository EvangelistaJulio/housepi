// individual.hpp //

#ifndef DEF_INDIVIDUAL__H
#define DEF_INDIVIDUAL__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "parameter.hpp"

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Community_Infec;
extern std::vector<double> Inst_Community_Infec;

/*
// TO BE UPDATED
//----------------------------------------------------------------------
// Augmented_Variable : Class declaration
//----------------------------------------------------------------------

class Augmented_Variable
{
	// Attributes
	double a_value = 0.0;
	double a_impossible = INFINITY;

public:
	// Getter
	double get() const { return a_value; }

	// Setter
	void init(std::mt19937_64 &rand_gen);
	void set(double new_value) { a_value = new_value; }

	// Computations
	double propose(std::mt19937_64 &rand_gen);
};
*/

//----------------------------------------------------------------------
// Individual : Class declaration
//----------------------------------------------------------------------

class Individual
{
	double i_augm_incub_period = 0.0;
	double i_augm_onset_time = INFINITY;

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
	double get_augm_incub_period() const { return i_augm_incub_period; }
	double get_augm_infec_time() const { return i_augm_onset_time - i_augm_incub_period; }
	double get_augm_onset_time() const { return i_augm_onset_time; }

	// Setter
	void init_augm_data(std::mt19937_64 &rand_gen);
	void set_augm_incub_period(double new_incub_period) { i_augm_incub_period = new_incub_period; }
	void set_augm_onset_time(double new_onset_time) { i_augm_onset_time = new_onset_time; }

	// Incubation period
	double log_incub_period_pdf(double x);

	// Propose augmented
	double propose_augm_incub_period(std::mt19937_64 &rand_gen);
	double propose_augm_onset_time(std::mt19937_64 &rand_gen);

	// Infectivity Profiles
	double infec_inst(double tinf, double study_period);
	double infec_cumul(double tinf, double study_period);
};

#endif //DEF_INDIVIDUAL__H