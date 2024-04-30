// individual.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

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
#define SUSCEPTIBLE 0
#define INFECTIOUS 1

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
	void assert_data(std::string dataID);

	//=====Rate Methods (e.g. Infectivity, Susceptibility, ...)=====//
	void assign_rate_per_group(std::string rateID, std::string parID);
	void _set_rate(std::string rateID, const double *x);
	double model_rate(std::string rateID) const;

	//=====Augmented Variables=====//
	void add_augm(std::string augmID, Augmented new_augm);
	void init_augm(std::string augmID, std::mt19937_64 &rng);
	void set_augm(std::string augmID, double new_value) { augm_var.at(augmID).set_value(new_value); }
	bool get_augm_infer(std::string augmID) const { return augm_var.at(augmID).get_infer(); }
	double get_augm(std::string augmID) const { return augm_var.at(augmID).get_value(); }

	//=====Augmented Variables: Model distribution=====//
	double augm_model_rand(std::string augmID, std::mt19937_64 &rng) { return augm_var.at(augmID).model_rand(rng); }
	double augm_model_log_pdf(std::string augmID) { return augm_var.at(augmID).model_log_pdf(augm_var.at(augmID).get_value()); }
	double augm_model_logratio(std::string augmID, double old_value, double new_value);

	//=====Augmented Variables: Proposer=====//
	double augm_propose(std::string augmID, std::mt19937_64 &rng);
	double augm_propose_logratio(std::string augmID, double old_value, double new_value);
	bool augm_valid(std::string augmID, double new_value) { return augm_var.at(augmID).valid(new_value); }
};

/////////////////////////////////////////////////////////
////////// Household : Class declaration ///////////////
////////////////////////////////////////////////////////

class Household
{
private:
	int house_index = -1;
	sdargs house_data = {};
	Individual **indiv_ptr = nullptr; // pointer to individauls' list

	// log-likelihood
	double loglike = -INF_DOUBLE;
	double buffer = -INF_DOUBLE;

	// accpetance rates for augmented variables
	std::unordered_map<std::string, int> num_accept_augm_moves = {};
	std::unordered_map<std::string, int> num_propos_augm_moves = {};

public:
	//=====Attributes=====//
	int house_size = 0;
	int first_case = -1;

	// list of house members for given status
	std::vector<int> infectious_list = {};
	std::vector<int> susceptibles_list = {};
	std::vector<int> unknown_status_list = {};

	// auxiliary variables to compute Bayesian diagnostics
	int dgn_ss = 0;
	long double dgn_sum_likelihood = 0.;
	long double dgn_sum_loglike = 0.;
	long double dgn_sum_sqr_loglike = 0.;

	//=====Constructors=====//
	Household() {}
	~Household(){};

	//======Setter======//
	void set_index(int index) { house_index = index; }
	void set_indiv_ptr(Individual **new_ptr) { indiv_ptr = new_ptr; }
	void link_individual(int i, int status);
	void change_status(int i, int new_status);
	void assign_data(std::string dataID, double value);

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

	void reset_augm_accept_rate(std::string augmID);
	void update_augm_accept_rate(std::string augmID, bool is_accepted);
	int get_augm_accept_counts(std::string augmID) { return num_accept_augm_moves[augmID]; }
	int get_augm_propos_counts(std::string augmID) { return num_propos_augm_moves[augmID]; }

	// update the values for the Bayesian diagnostics
	void update_diagnostics();
};

#endif // INDIVIDUAL__H