// individual.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#ifndef INDIVIDUAL__H
#define INDIVIDUAL__H

#include <string>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <cmath>

#include "varia.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Individual : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Individual
{
private:
	//=====Attributes (private)=====//
	int status = 0;	  // status of the individual
	int index = 0;	  // database index of the individual
	sdargs data = {}; // database associated to individual

public:
	//=====Attributes (public)=====//
	std::unordered_map<std::string, std::string> which_rate = {};

	//=====Constructors=====//
	Individual(sdargs data);
	Individual() = default;
	~Individual(){};

	//=====Getter=====//
	int get_status() const { return status; }
	void set_status(int d) { status = d; }
	int get_index() const { return index; }

	// Note : get_data is bottleneck if called every iteration.
	// Deprecated outisde of pre-processing.
	double get_data(std::string dataID) const;
	void assert_data(std::string dataID);

	//=====Rate Methods (e.g. Infectivity, Susceptibility, ...)=====//
	void assign_rate_per_group(std::string rateID, std::string parID);
	void set_rate(std::string rateID, const double *x) { rate_ptr.emplace(rateID, x); }
	double get_group_rate(std::string rateID) const { return *rate_ptr.at(rateID); }

	//=====Augmented Variables=====//
	void assign_augm(std::string augmID, Augmented new_augm);
	void init_augm(std::string augmID, sys_rnd_gen &rng);
	void set_augm(std::string augmID, double new_value) { augm_var.at(augmID).set_value(new_value); }
	bool get_augm_infer(std::string augmID) const { return augm_var.at(augmID).get_infer(); }
	double get_augm(std::string augmID) const { return augm_var.at(augmID).get_value(); }

	//-----Augmented Variables: Model distribution-----//
	double augm_model_rand(std::string augmID, sys_rnd_gen &rng) { return augm_var.at(augmID).model_rand(rng); }
	double augm_model_log_pdf(std::string augmID) { return augm_var.at(augmID).model_log_pdf(augm_var.at(augmID).get_value()); }
	double augm_model_logratio(std::string augmID, double old_value, double new_value);

	//-----Augmented Variables: Proposer-----//
	double augm_propose(std::string augmID, sys_rnd_gen &rng);
	double augm_propose_logratio(std::string augmID, double old_value, double new_value);
	bool augm_valid(std::string augmID, double new_value) { return augm_var.at(augmID).valid(new_value); }

private:
	sdargs_ext rate_ptr = {};
	std::unordered_map<std::string, Augmented> augm_var = {};
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Household
{
private:
	//=====Attributes (private)=====//
	int house_index = -1;		  // database index of the household
	sdargs house_data = {};		  // database associated to individual

	double loglike = -INF_DOUBLE; // log-likelihood value
	double buffer = -INF_DOUBLE;  // store buffer value of the log-likelihood

	// accpetance rates for augmented variables
	std::unordered_map<std::string, int> num_accept_augm_moves = {};
	std::unordered_map<std::string, int> num_propos_augm_moves = {};

public:
	//=====Attributes (public)=====//
	std::vector<Individual> indiv = {};
	int house_size = 0;
	int num_infectious = 0;
	int num_susceptible = 0;
	int num_missing_status = 0;
	int first_case = -1;

	// contact rate address
	std::vector<std::vector<const double *>> contact_rates;
	std::vector<std::vector<std::string>> which_rate;

	// wheter members correspond to a given status or not
	std::vector<bool> is_infectious = {};
	std::vector<bool> is_susceptible = {};
	std::vector<bool> is_missing_status = {};

	// list with ONLY indexes of members with a given status
	std::set<int> infectious_set = {};
	std::set<int> susceptible_set = {};
	std::set<int> missing_status_set = {};

	//===== Diagnostics =====//
	int dgn_sample_size = 0;		  // sample size
	long double dgn_likelihood = 0.;  // mean (house) likelihood
	long double dgn_loglike = 0.;	  // mean (house) log-likelihood
	long double dgn_sqr_loglike = 0.; // mean squared (house) log-likelihood

	//===== Constructors =====//
	Household() {}
	~Household(){};

	//====== Setter ======//
	void assign_individual(Individual this_indiv);
	void set_index(int index) { house_index = index; }
	void assign_data(std::string dataID, double value);
	void assign_rate_per_contact(int from, int to, std::string parID) { which_rate[from][to] = parID; }
	void change_status(int i, int new_status);

	//====== Getter ======//
	int get_index() const { return house_index; }
	double get_data(std::string dataID) const;
	double get_group_rate(int i, std::string rateID) const { return indiv[i].get_group_rate(rateID); }

	//====== Contact Rates ======//
	void set_rate(int from, int to, const double *x) { contact_rates[from][to] = x; }
	double get_contact_rate(int from, int to) const { return (*contact_rates[from][to]); }

	//====== Loglikelihood ======//
	double get_loglike() const { return loglike; }
	void set_loglike(double new_loglike)
	{
		buffer = loglike;
		loglike = new_loglike;
	}
	void refuse() { loglike = buffer; }
	void init_diagnostics(int sample_size);
	void update_diagnostics();

	//===== Acceptation Rate =====//
	void reset_augm_accept_rate(std::string augmID);
	void update_augm_accept_rate(std::string augmID, bool is_accepted);
	int get_augm_accept_counts(std::string augmID) { return num_accept_augm_moves[augmID]; }
	int get_augm_propos_counts(std::string augmID) { return num_propos_augm_moves[augmID]; }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#endif // INDIVIDUAL__H