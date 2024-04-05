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
	int home = -1;
	std::unordered_map<std::string, std::string> _which_rate = {};

	//=====Constructors=====//
	Individual(sdargs data);
	Individual() = default;
	~Individual(){};

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
	int house_index = -1;
	Individual **indiv_ptr = nullptr;
	sdargs house_data = {};
	double loglike = -INFINITY;
	double buffer = -INFINITY;

	// auxiliary variables to compute Bayesian diagnostics
	int dgn_sample_size = 0;
	double dgn_mean_ll = 0.;
	double dgn_mean_sqr_ll = 0.;

public:
	//=====Attributes=====//
	int house_size = 0;
	int first_case = -1;
	std::vector<int> infectious_list = {};
	std::vector<int> susceptibles_list = {};
	std::vector<int> unknown_status_list = {};

	//=====Constructors=====//
	Household() {}
	~Household(){};

	//======Setter======//
	void set_index(int index) { house_index = index; }
	void set_indiv_ptr(Individual **new_ptr) { indiv_ptr = new_ptr; }
	void link_individual(int i, int status);
	void change_status(int i, int new_status);
	void assign_data(std::string dataID, double value);
	void set_data(std::string dataID, double value) { house_data.at(dataID) = value; }

	//======Getter======//
	int get_index() const { return house_index; }
	std::vector<int> get_house_members();

	/* get_data : bottleneck if data are called every iteration. Deprecated outisde of pre-processing.*/
	double get_data(std::string dataID) const;

	//======Loglike======//
	double get_loglike() const { return loglike; }
	void set_loglike(double new_loglike)
	{
		buffer = loglike;
		loglike = new_loglike;
	}
	void refuse() { loglike = buffer; }

	// update the values for the Bayesian diagnostics
	void update_dgn() 
	{
		dgn_mean_ll = (loglike + dgn_mean_ll * dgn_sample_size) / (dgn_sample_size + 1);
		dgn_mean_sqr_ll = (std::pow(loglike, 2) + dgn_mean_sqr_ll * dgn_sample_size) / (dgn_sample_size + 1);
		dgn_sample_size ++;
	}
};

#endif // INDIVIDUAL__H