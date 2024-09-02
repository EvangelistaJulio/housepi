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

	//===== Constructors =====//
	Individual() = default;
	Individual(sdargs data);
	~Individual() {};

	//===== Basic Methods =====//
	int get_status() const { return status; }
	int get_index() const { return index; }

	void set_status(int s) { status = s; }
	void set_index(int i) { index = i; }

	void disable_assert();
	void assert_data(const std::string d) const;
	double get_data(const std::string d) const { return (this->*_get_data_method_)(d); }

	//=====Rate Methods (e.g. Infectivity, Susceptibility, ...)=====//
	void assign_rate_per_group(const std::string r, const std::string parID);
	void set_rate(const std::string r, const double *x) { rate_ptr.emplace(r, x); }
	double get_group_rate(const std::string r) const { return *rate_ptr.at(r); }

	//=====Augmented Variables=====//
	void assert_augm(const std::string a) const;
	void assign_augm(const std::string a, double value,
					 const std::string model_name, sdargs_ext model_args_ext,
					 const std::string rule_name, sdargs rule_args,
					 double rule_from, sdargs value_range);
	void init_augm(const std::string a, sys_rnd_gen &prng);
	void set_augm(const std::string a, double new_value) { augm_var.at(a).set_value(new_value); }
	bool get_augm_infer(const std::string a) const { return augm_var.at(a).get_infer(); }
	double get_augm(const std::string a) const { return (this->*_get_augm_method_)(a); }

	//-----Augmented Variables: Model distribution-----//
	double augm_model_rand(const std::string a, sys_rnd_gen &prng) { return augm_var.at(a).model_rand(prng); }
	double augm_model_log_pdf(const std::string a) { return augm_var.at(a).model_log_pdf(augm_var.at(a).get_value()); }
	double augm_model_logratio(const std::string a, double old_value, double new_value);

	//-----Augmented Variables: Proposer-----//
	double augm_propose(const std::string a, sys_rnd_gen &prng);
	double augm_propose_logratio(const std::string a, double old_value, double new_value);
	bool augm_valid(const std::string a, double new_value) { return augm_var.at(a).valid(new_value); }

private:
	sdargs_ext rate_ptr = {};
	std::unordered_map<std::string, Augmented> augm_var = {};

	//===== Setter (auxiliary methods) =====//
	void (Individual::*_set_augm_method_)(const std::string, double) = &Individual::_set_augm_assert_;
	void _set_augm_assert_(const std::string a, double new_value);
	void _set_augm_default_(const std::string a, double new_value);

	//===== Getter (auxiliary methods) =====//
	double (Individual::*_get_augm_method_)(const std::string) const = &Individual::_get_augm_assert_;
	double _get_augm_assert_(const std::string a) const;
	double _get_augm_default_(const std::string a) const;

	double (Individual::*_get_data_method_)(const std::string) const = &Individual::_get_data_assert_;
	double _get_data_assert_(const std::string d) const;
	double _get_data_default_(const std::string d) const;
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
	int house_index = -1;	// database index of household
	sdargs house_data = {}; // database associated to household (manual)

	double house_loglike = -INF_DOUBLE; // current loglikelihood value
	double buffer = -INF_DOUBLE;		// buffer loglikelihood value

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
	ulong dgn_ss = 0;		   // sample size
	long double dgn_like = 0.; // average (house) likelihood
	long double dgn_avg = 0.;  // average (house) loglikelihood
	long double dgn_sqr = 0.;  // average squared (house) loglikelihood
	long double dgn_var() { return (dgn_sqr - std::pow(dgn_avg, 2)) * dgn_ss / (dgn_ss - 1); }

	//===== Constructors =====//
	Household() {}
	~Household() {};

	//====== Setter ======//
	void assign_individual(Individual this_indiv);
	void set_index(int index) { house_index = index; }
	void assign_data(const std::string d, double value);
	void assign_rate_per_contact(int from, int to, const std::string parID) { which_rate[from][to] = parID; }
	void change_status(int i, int new_status);

	//====== Getter ======//
	int get_index() const { return house_index; }
	void assert_data(const std::string d) const;
	double get_data(const std::string d) const;
	double get_group_rate(int i, const std::string r) const { return indiv[i].get_group_rate(r); }

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
	void init_diagnostics(int sample_size);
	void update_diagnostics();

	//===== Acceptation Rate =====//
	void reset_augm_accept_rate(const std::string a);
	void update_augm_accept_rate(const std::string a, bool is_accepted);
	int get_augm_accept_counts(const std::string a) { return num_accept_augm_moves[a]; }
	int get_augm_propos_counts(const std::string a) { return num_propos_augm_moves[a]; }
};

#endif // INDIVIDUAL__H