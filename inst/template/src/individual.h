/*
 * individual.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef INDIVIDUAL__H
#define INDIVIDUAL__H

#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

#include "varia.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// status aliases
constexpr int INFECTIOUS = 1;	   // status of infectious (>=)
constexpr int SUSCEPTIBLE = 0;	   // status of susceptible (==)
constexpr int MISSING_STATUS = -1; // missing status (==)

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Augmented : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Augmented : public Base_Variable
{
public:
	//===== Constructor =====//
	Augmented(const std::string &augm_name = UNDEFINED, double value = std::numeric_limits<double>::quiet_NaN(),
			  const std::string &model_name = UNDEFINED, dictargs_ptr model_args_ptr = {}, dictargs range = {},
			  const std::string &rule_name = UNDEFINED, dictargs rule_args = {}, double rule_from = INF_DOUBLE) : Base_Variable(augm_name, value)
	{
		assign_distribution(model_name, model_args_ptr, range);
		assign_proposer(rule_name, rule_args, rule_from);
	}
};

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
	std::string ID = UNDEFINED; // user ID
	int status = 0;				// status of the individual
	dictargs data = {};			// database associated to individual

	// NOTE : for the moment, I can't retrieve house from individual

public:
	//=====Attributes (public)=====//
	std::unordered_map<std::string, std::string> which_rate = {};

	//===== Constructors =====//
	Individual() = default;
	Individual(const std::string &new_ID, int new_status, const dictargs &data);
	~Individual() {};

	//===== Getter Methods =====//
	std::string get_ID() const { return ID; }
	int get_status() const { return status; }
	double get_data(const std::string &dataID) const { return (this->*_get_data_method_)(dataID); }
	const double *get_data_pointer(const std::string &dataID) const { return (this->*_get_data_pointer_method_)(dataID); }

	//===== Assert Methods =====//
	void disable_assert();
	void assert_data(const std::string &dataID) const;

	//===== Setter Methods =====//
	void set_ID(const std::string &new_ID) { ID = new_ID; }
	void set_status(int s) { status = s; }

	//===== Rate Methods (e.g. Infectivity, Susceptibility, ...) =====//
	void assign_group_rate(const std::string &rateID, const std::string &parID);
	void set_rate(const std::string &rateID, const double *x) { rate_ptr.emplace(rateID, x); }
	double get_group_rate(const std::string &rateID) const { return *rate_ptr.at(rateID); }

	//===== Augmented Variables =====//
	void assign_augm(const Augmented &new_augm) { augm_var[new_augm.get_name()] = new_augm; }
	void assert_augm(const std::string &augmID) const;
	void augm_model_rand(const std::string &augmID, sys_rnd_gen &prng) { augm_var.at(augmID).rand(prng); }
	void set_augm(const std::string &augmID, double new_value) { augm_var.at(augmID).set_value(new_value); }
	bool get_augm_infer(const std::string &augmID) const { return augm_var.at(augmID).get_infer(); }
	double get_augm(const std::string &augmID) const { return (this->*_get_augm_method_)(augmID); }

	//===== Metropolis-Hastings =====//
	void augm_propose(const std::string &augmID, sys_rnd_gen &prng) { augm_var.at(augmID).propose(prng); }
	bool augm_valid_move(const std::string &augmID) const { return augm_var.at(augmID).valid_move(); }
	bool propose_valid_augm_move(const std::string &augmID, sys_rnd_gen &prng);

	double augm_model_logpdf(const std::string &augmID) const { return augm_var.at(augmID).distr_logpdf(); }
	double augm_model_logratio(const std::string &augmID) const { return augm_var.at(augmID).distr_logratio(); }
	double augm_propose_logratio(const std::string &augmID) const { return augm_var.at(augmID).propose_logratio(); }
	void augm_restore_value(const std::string &augmID) { augm_var.at(augmID).restore_value(); }
	double augm_move_logratio(const std::string &augmID) const { return augm_var.at(augmID).get_move_logratio(); };

private:
	dictargs_ptr rate_ptr = {};
	std::unordered_map<std::string, Augmented> augm_var = {};

	//===== Setter (auxiliary methods) =====//
	void (Individual::*_set_augm_method_)(const std::string &, double) = &Individual::_set_augm_assert_;
	void _set_augm_assert_(const std::string &augmID, double new_value);
	void _set_augm_default_(const std::string &augmID, double new_value);

	//===== Getter (auxiliary methods) =====//
	double (Individual::*_get_augm_method_)(const std::string &) const = &Individual::_get_augm_assert_;
	double _get_augm_assert_(const std::string &augmID) const;
	double _get_augm_default_(const std::string &augmID) const;

	double (Individual::*_get_data_method_)(const std::string &) const = &Individual::_get_data_assert_;
	double _get_data_assert_(const std::string &dataID) const;
	double _get_data_default_(const std::string &dataID) const;

	const double *(Individual::*_get_data_pointer_method_)(const std::string &) const = &Individual::_get_data_pointer_assert_;
	const double *_get_data_pointer_assert_(const std::string &dataID) const;
	const double *_get_data_pointer_default_(const std::string &dataID) const;
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
	std::string ID = UNDEFINED; // user ID
	int index = -1;				// index of household
	int house_size = 0;			// number of house members

	double house_loglike = -INF_DOUBLE; // current loglikelihood value
	double buffer = -INF_DOUBLE;		// buffer loglikelihood value

	// accpetance rates for augmented variables
	std::unordered_map<std::string, int> num_accept_augm_moves = {};
	std::unordered_map<std::string, int> num_propos_augm_moves = {};

	// wheter members correspond to a given status or not
	std::vector<bool> _is_infectious_ = {};
	std::vector<bool> _is_susceptible_ = {};
	std::vector<bool> _is_missing_status_ = {};

public:
	//=====Attributes (public)=====//
	std::vector<Individual> indiv = {};
	int num_infectious = 0;
	int num_susceptible = 0;
	int num_missing_status = 0;

	// contact rate pointers
	std::vector<std::vector<const double *>> contact_rates;
	std::vector<std::vector<std::string>> which_rate;

	// list with ONLY indexes of members with a given status
	std::set<int> infectious_set = {};
	std::set<int> susceptible_set = {};
	std::set<int> missing_status_set = {};

	//===== Constructors =====//
	Household() {}
	~Household() {};

	//====== Setter ======//
	void assign_individual(Individual this_indiv);
	void set_ID(const std::string &new_ID) { ID = new_ID; }
	void set_index(int h) { index = h; }
	void assign_rate_per_contact(int from, int to, const std::string &parID) { which_rate[from][to] = parID; }
	bool change_status(int i, int new_status);

	//====== Getter ======//
	std::string get_ID() const { return ID; }
	int get_index() const { return index; }
	int size() const { return house_size; }

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
	std::vector<double> dgn_loglike = {};
	void init_diagnostics(ulong sample_size);
	void update_diagnostics();
	std::vector<double> get_dgn_loglike() const { return dgn_loglike; }

	//===== Acceptation Rate =====//
	void reset_augm_accept_rate(const std::string &augmID);
	void update_augm_accept_rate(const std::string &augmID, bool is_accepted);
	int get_augm_accept_counts(const std::string &augmID) { return num_accept_augm_moves[augmID]; }
	int get_augm_propos_counts(const std::string &augmID) { return num_propos_augm_moves[augmID]; }
};

#endif // INDIVIDUAL__H