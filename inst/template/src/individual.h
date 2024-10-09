/*
 * individual.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef INDIVIDUAL__H
#define INDIVIDUAL__H

#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

#include "varia.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// status aliases
constexpr int INFECTIOUS = 1;      // status of infectious (>=)
constexpr int SUSCEPTIBLE = 0;     // status of susceptible (==)
constexpr int MISSING_STATUS = -1; // missing status (==)

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
	void set_status(int s) { status = s; }
	void disable_assert();
	void assert_data(const string &dataID) const;
	double get_data(const string &dataID) const { return (this->*_get_data_method_)(dataID); }
	const double *get_data_pointer(const string &dataID) const { return (this->*_get_data_pointer_method_)(dataID); }

	//=====Rate Methods (e.g. Infectivity, Susceptibility, ...)=====//
	void assign_group_rate(const string &rateID, const string &parID);
	void set_rate(const string &rateID, const double *x) { rate_ptr.emplace(rateID, x); }
	double get_group_rate(const string &rateID) const { return *rate_ptr.at(rateID); }

	//===== Augmented Variables =====//
	void assign_augm(const Augmented &new_augm) { augm_var[new_augm.get_name()] = new_augm; }
	void assert_augm(const string &augmID) const;
	void augm_model_rand(const string &augmID, sys_rnd_gen &prng) { augm_var.at(augmID).rand(prng); }
	void set_augm(const string &augmID, double new_value) { augm_var.at(augmID).set_value(new_value); }
	bool get_augm_infer(const string &augmID) const { return augm_var.at(augmID).get_infer(); }
	double get_augm(const string &augmID) const { return (this->*_get_augm_method_)(augmID); }

	//===== Metropolis-Hastings =====//
	void augm_propose(const string &augmID, sys_rnd_gen &prng) { augm_var.at(augmID).propose(prng); }
	bool augm_valid_move(const string &augmID) const { return augm_var.at(augmID).valid_move(); }
	bool propose_valid_augm_move(const string &augmID, sys_rnd_gen &prng);

	double augm_model_logpdf(const string &augmID) const { return augm_var.at(augmID).distr_logpdf(); }
	double augm_model_logratio(const string &augmID) const { return augm_var.at(augmID).distr_logratio(); }
	double augm_propose_logratio(const string &augmID) const { return augm_var.at(augmID).propose_logratio(); }
	void augm_restore_value(const string &augmID) { augm_var.at(augmID).restore_value(); }
	double augm_move_logratio(const string &augmID) const { return augm_var.at(augmID).get_move_logratio(); };

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
	//===== Attributes (private) =====//
	int index = -1;			  // index of household
	dictargs house_data = {}; // database associated to household (manual)
	int house_size = 0;

	double house_loglike = -INF_DOUBLE; // current loglikelihood value
	double buffer = -INF_DOUBLE;		// buffer loglikelihood value

	// accpetance rates for augmented variables
	std::unordered_map<string, int> num_accept_augm_moves = {};
	std::unordered_map<string, int> num_propos_augm_moves = {};

	// wheter members correspond to a given status or not
	vector<bool> _is_infectious_ = {};
	vector<bool> _is_susceptible_ = {};
	vector<bool> _is_missing_status_ = {};

public:
	//=====Attributes (public)=====//
	vector<Individual> indiv = {};
	int num_infectious = 0;
	int num_susceptible = 0;
	int num_missing_status = 0;

	// contact rate pointers
	vector<vector<const double *>> contact_rates;
	vector<vector<string>> which_rate;

	// list with ONLY indexes of members with a given status
	std::set<int> infectious_set = {};
	std::set<int> susceptible_set = {};
	std::set<int> missing_status_set = {};

	//===== Constructors =====//
	Household() {}
	~Household() {};

	//====== Setter ======//
	void assign_individual(Individual this_indiv);
	void set_id(int id) { index = id; }
	void assign_data(const string &dataID, double value);
	void assign_rate_per_contact(int from, int to, const string &parID) { which_rate[from][to] = parID; }
	bool change_status(int i, int new_status);

	//====== Getter ======//
	int size() const { return house_size; }
	int get_id() const { return index; }
	void assert_data(const string &dataID) const;
	double get_data(const string &dataID) const;

	// wheter members correspond to a given status or not
	bool is_infectious(int i) const { return _is_infectious_[i]; }
	bool is_susceptible(int i) const { return _is_susceptible_[i]; }
	bool is_missing_status(int i) const { return _is_missing_status_[i]; }

	//====== Contact Rates ======//
	void set_contact_rate(int from, int to, const double *x) { contact_rates[from][to] = x; }
	double get_contact_rate(int from, int to) const { return (*contact_rates[from][to]); }

	//====== Loglikelihood ======//
	double get_loglike() const { return house_loglike; }
	void set_loglike(double new_loglike)
	{
		buffer = house_loglike;
		house_loglike = new_loglike;
	}
	void restore_loglike() { house_loglike = buffer; }
	double get_delta_loglike() const { return house_loglike - buffer; }

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