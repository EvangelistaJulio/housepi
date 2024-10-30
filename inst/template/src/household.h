/*
 * household.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEHOLD__H
#define HOUSEHOLD__H

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

#include "individual.h"

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
public:
	//===== ATTRIBUTES =====//
	std::vector<Individual> indiv = {};
	
	//===== CONSTRUCTOR =====//
	Household();
	Household(const std::string &ID, int index);
	~Household();

	//===== SETTER =====//
	void assign_individual(const std::string &ID, int status, const dictargs &data);
	bool change_status(int i, int new_status);
	void init_contact_matrix();
	void assign_rate_per_contact(int from, int to, const std::string &parID);
	void set_contact_rate(int from, int to, const double *x);

	//===== GETTER =====//
	const std::string &get_ID() const;
	int get_index() const;
	int size() const;
	const std::string &get_group_rate_map(int from, int to) const;
	double get_contact_rate(int from, int to) const;
	bool is_infectious(int i) const;
	bool is_susceptible(int i) const;
	bool is_missing_status(int i) const;
	int get_num_infectious() const;
	int get_num_susceptible() const;
	int get_num_missing_status() const;
	const std::set<int> &get_infectious_set() const;
	const std::set<int> &get_susceptible_set() const;
	const std::set<int> &get_missing_status_set() const;
	
	//===== LOGLIKELIHOOD =====//

	//----- Setter -----//
	void set_loglike(double x);
	void restore_loglike();

	//----- Getter -----//
	double get_loglike() const;
	double get_delta_loglike() const;

	//===== DIAGNOSTICS =====//
	void init_diagnostics(ulong sample_size);
	void update_diagnostics();
	std::vector<double> get_dgn_loglike() const { return dgn_loglike; }

private:
	//===== Attributes =====//
	std::string ID = UNDEFINED; // user ID
	int index = -1;				// index of household
	int house_size = 0;			// number of house members

	double house_loglike = -INF_DOUBLE; // current loglikelihood value
	double buffer = -INF_DOUBLE;		// buffer loglikelihood value

	int num_infectious = 0;
	int num_susceptible = 0;
	int num_missing_status = 0;

	// wheter members correspond to a given status or not
	std::vector<bool> _is_infectious_ = {};
	std::vector<bool> _is_susceptible_ = {};
	std::vector<bool> _is_missing_status_ = {};

	ulong dgn_count = 0;
	std::vector<double> dgn_loglike = {};

	// contact rate pointers
	std::vector<std::vector<const double *>> contact_rates;
	std::vector<std::vector<std::string>> _which_rate_;

	// list with ONLY indexes of members with a given status
	std::set<int> infectious_set = {};
	std::set<int> susceptible_set = {};
	std::set<int> missing_status_set = {};
};

#endif // HOUSEHOLD__H