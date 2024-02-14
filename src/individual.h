// individual.h //

#ifndef INDIVIDUAL__H
#define INDIVIDUAL__H

#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "varia.h"

// model
#define HEALTHY 0
#define SYMPTOMATIC 1

/////////////////////////////////////////////////////////
///////// Individual : Class declaration ///////////////
////////////////////////////////////////////////////////

class Individual
{
private:
	//=====Private Attributes=====//
	int status = 0;
	int index = 0;

	sdargs data = {};
	std::unordered_map<std::string, Augmented> augm_var = {};
	sdargs_ext _rate_ptr = {};

public:
	//=====Public Attributes=====//
	size_t home = -1;
	std::unordered_map<std::string, std::string> _which_rate = {};

	//=====Constructors=====//
	Individual(sdargs data);
	Individual() = default;
	~Individual(){};

	//=====Methods for the Model=====//
	double get_infec_time() const;

	//=====Getter=====//
	bool get_status() const { return status; }
	void set_status(int d) { status = d; }
	int get_index() const { return index; }

	/* get_data : bottleneck if data are called every iteration. Deprecated outisde of pre-processing.*/
	double get_data(std::string dataID) const;

	//=====Rate Methods (e.g. Infectivity, Susceptibility, ...)=====//
	void assign_rate_per_group(std::string rateID, std::string parID);
	void _set_rate(std::string rateID, const double *x);
	double model_rate(std::string rateID) const;

	//=====Augmented Variables=====//
	void add_augm(std::string augmID, Augmented new_augm);
	void add_augm(std::string augmID, Augmented new_augm, std::mt19937_64 &rng);
	double get_augm(std::string augmID) const { return augm_var.at(augmID).get_value(); }
	bool get_augm_infer(std::string augmID) { return augm_var.at(augmID).get_infer(); }
	void set_augm(std::string augmID, double new_val) { augm_var.at(augmID).set_value(new_val); }
	void init_augm(std::string augmID, std::mt19937_64 &rng);

	//=====Augmented Variables: Model distribution=====//
	double augm_model_rand(std::string augmID, std::mt19937_64 &rng) { return augm_var.at(augmID).model_rand(rng); }
	double augm_model_log_pdf(std::string augmID);
	double augm_model_log_ratio(std::string augmID, double old_value, double new_val);

	//=====Augmented Variables: Proposer=====//
	double augm_propose(std::string augmID, std::mt19937_64 &rng);
	double augm_propose_log_ratio(std::string augmID, double old_value, double new_val);
	bool augm_valid(std::string augmID, double new_val) { return augm_var.at(augmID).valid(new_val); }
};

/////////////////////////////////////////////////////////
////////// Household : Class declaration ///////////////
////////////////////////////////////////////////////////

class Household
{
private:
	Individual **indiv_ptr = nullptr;
	double loglike = -INFINITY;
	double buffer = -INFINITY;

public:
	//=====Attributes=====//
	int house_index = -1;
	size_t size = 0;
	size_t first_case = -1;
	double start_followup = INFINITY;
	double end_followup = 0;
	std::vector<size_t> susceptibles_list = {};
	std::vector<size_t> infectors_list = {};

	//=====Constructors=====//
	Household() {}
	~Household(){};

	//======Getter======//
	int get_index() const { return house_index; }

	//======Setter======//
	void set_index(int index) { house_index = index; }
	void set_indiv_ptr(Individual **new_ptr) { indiv_ptr = new_ptr; }
	void link_individual(size_t i, int status);
	void change_status(size_t i, int new_status);
	void update_first_case();
	void assign_end_followup();

	//======Loglike======//
	double get_loglike() const { return loglike; }
	void set_loglike(double new_loglike)
	{
		buffer = loglike;
		loglike = new_loglike;
	}
	void refuse() { loglike = buffer; }
};

#endif // INDIVIDUAL__H