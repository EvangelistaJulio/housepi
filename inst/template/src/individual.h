/*
 * individual.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef INDIVIDUAL__H
#define INDIVIDUAL__H

#include <algorithm>
#include <set>
#include <unordered_map>
#include <cmath>

#include "varia.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

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
	int status = 0;		// status of the individual
	int index = 0;		// database index of the individual
	dictargs data = {}; // database associated to individual

public:
	//=====Attributes (public)=====//
	std::unordered_map<string, string> which_rate = {};

	//===== Constructors =====//
	Individual() = default;
	Individual(dictargs data);
	~Individual() {};

	//===== Basic Methods =====//
	int get_status() const { return status; }
	int get_index() const { return index; }

	void set_status(int s) { status = s; }
	void set_index(int i) { index = i; }

	void disable_assert();
	void assert_data(const string &dataID) const;
	double get_data(const string &dataID) const { return (this->*_get_data_method_)(dataID); }
	const double *get_data_pointer(const string &dataID) const { return (this->*_get_data_pointer_method_)(dataID); }

	//=====Rate Methods (e.g. Infectivity, Susceptibility, ...)=====//
	void assign_rate_per_group(const string &rateID, const string &parID);
	void set_rate(const string &rateID, const double *x) { rate_ptr.emplace(rateID, x); }
	double get_group_rate(const string &rateID) const { return *rate_ptr.at(rateID); }

	//=====Augmented Variables=====//
	void assign_augm(const string &augmID, double value,
					 const string &model_name, const dictargs_ptr &model_args,
					 const string &rule_name, const dictargs &rule_args,
					 const dictargs &range, double rule_from);
	void assert_augm(const string &augmID) const;
	void init_augm(const string &augmID, sys_rnd_gen &prng);
	void set_augm(const string &augmID, double new_value) { augm_var.at(augmID).set_value(new_value); }
	bool get_augm_infer(const string &augmID) const { return augm_var.at(augmID).get_infer(); }
	double get_augm(const string &augmID) const { return (this->*_get_augm_method_)(augmID); }

	//-----Augmented Variables: Model distribution-----//
	double augm_model_rand(const string &augmID, sys_rnd_gen &prng) { return augm_var.at(augmID).model_rand(prng); }
	double augm_model_log_pdf(const string &augmID) const { return augm_var.at(augmID).model_log_pdf(augm_var.at(augmID).get_value()); }
	double augm_model_logratio(const string &augmID, double old_value, double new_value) const;

	//-----Augmented Variables: Proposer-----//
	double augm_propose(const string &augmID, sys_rnd_gen &prng);
	double augm_propose_logratio(const string &augmID, double old_value, double new_value) const;
	bool augm_valid(const string &augmID, double new_value) const { return augm_var.at(augmID).valid(new_value); }

private:
	dictargs_ptr rate_ptr = {};
	std::unordered_map<string, Augmented> augm_var = {};

	//===== Setter (auxiliary methods) =====//
	void (Individual::*_set_augm_method_)(const string &, double) = &Individual::_set_augm_assert_;
	void _set_augm_assert_(const string &augmID, double new_value);
	void _set_augm_default_(const string &augmID, double new_value);

	//===== Getter (auxiliary methods) =====//
	double (Individual::*_get_augm_method_)(const string &) const = &Individual::_get_augm_assert_;
	double _get_augm_assert_(const string &augmID) const;
	double _get_augm_default_(const string &augmID) const;

	double (Individual::*_get_data_method_)(const string &) const = &Individual::_get_data_assert_;
	double _get_data_assert_(const string &dataID) const;
	double _get_data_default_(const string &dataID) const;

	const double *(Individual::*_get_data_pointer_method_)(const string &) const = &Individual::_get_data_pointer_assert_;
	const double *_get_data_pointer_assert_(const string &dataID) const;
	const double *_get_data_pointer_default_(const string &dataID) const;
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
	int house_index = -1;	  // database index of household
	dictargs house_data = {}; // database associated to household (manual)

	double house_loglike = -INF_DOUBLE; // current loglikelihood value
	double buffer = -INF_DOUBLE;		// buffer loglikelihood value

	// accpetance rates for augmented variables
	std::unordered_map<string, int> num_accept_augm_moves = {};
	std::unordered_map<string, int> num_propos_augm_moves = {};

public:
	//=====Attributes (public)=====//
	vector<Individual> indiv = {};
	int house_size = 0;
	int num_infectious = 0;
	int num_susceptible = 0;
	int num_missing_status = 0;

	// contact rate pointers
	vector<vector<const double *>> contact_rates;
	vector<vector<string>> which_rate;

	// wheter members correspond to a given status or not
	vector<bool> is_infectious = {};
	vector<bool> is_susceptible = {};
	vector<bool> is_missing_status = {};

	// list with ONLY indexes of members with a given status
	std::set<int> infectious_set = {};
	std::set<int> susceptible_set = {};
	std::set<int> missing_status_set = {};

	//===== Constructors =====//
	Household() {}
	~Household() {};

	//====== Setter ======//
	void assign_individual(Individual this_indiv);
	void set_index(int index) { house_index = index; }
	void assign_data(const string &dataID, double value);
	void assign_rate_per_contact(int from, int to, const string &parID) { which_rate[from][to] = parID; }
	void change_status(int i, int new_status);

	//====== Getter ======//
	int get_index() const { return house_index; }
	void assert_data(const string &dataID) const;
	double get_data(const string &dataID) const;
	double get_group_rate(int i, const string &rateID) const { return indiv[i].get_group_rate(rateID); }

	//====== Contact Rates ======//
	void set_rate(int from, int to, const double *x) { contact_rates[from][to] = x; }
	double get_contact_rate(int from, int to) const { return (*contact_rates[from][to]); }

	//====== Loglikelihood ======//
	double get_loglike() const { return house_loglike; }
	void set_loglike(double new_loglike)
	{
		buffer = house_loglike;
		house_loglike = new_loglike;
	}
	void refuse() { house_loglike = buffer; }

	//===== Diagnostics =====//
	ulong dgn_count = 0;
	vector<double> dgn_loglike = {};
	void init_diagnostics(ulong sample_size);
	void update_diagnostics();
	vector<double> get_dgn_loglike() const { return dgn_loglike; }

	//===== Acceptation Rate =====//
	void reset_augm_accept_rate(const string &augmID);
	void update_augm_accept_rate(const string &augmID, bool is_accepted);
	int get_augm_accept_counts(const string &augmID) { return num_accept_augm_moves[augmID]; }
	int get_augm_propos_counts(const string &augmID) { return num_propos_augm_moves[augmID]; }
};

#endif // INDIVIDUAL__H