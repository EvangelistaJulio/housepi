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

#include "tools.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class Household_Impl;

class Household
{
public:
	//===== CONSTRUCTOR =====//
	Household();
	Household(const std::string &ID, int index);
	Household(Household&& other) noexcept;        // Move constructor
    Household& operator=(Household&& other) noexcept; // Move assignment
	~Household();

	//===== SETTER =====//
	void assign_individual(const std::string &ID, int status, const dictargs &data);
	bool change_status(int i, int new_status);
	void init_contact_matrix();
	void assign_contact_rate(int from, int to, const std::string &parID);
	void set_contact_rate(int from, int to, const double *x);
	void disable_assert();

	//===== GETTER =====//
	const std::string house_ID() const;
	int index() const;
	int size() const;
	int num_infectious() const;
	int num_susceptible() const;
	int num_missing_status() const;
	const std::set<int> &infectious_set() const;
	const std::set<int> &susceptible_set() const;
	const std::set<int> &missing_status_set() const;

	//===== INDIVIDUAL =====//

	//----- Individual : Setter -----//
	void assign_group_rate(int i, const std::string &rateID, const std::string &parID);
	void set_rate(int i, const std::string &rateID, const double *ptr);

	//----- Individual : Getter -----//
	const std::string ID(int i) const;
	double data(int i, const std::string &dataID) const;
	const double *data_pointer(int i, const std::string &dataID) const;
	std::unordered_map<std::string, std::string> get_group_rate_map(int i) const;
	double group_rate(int i, const std::string &rateID) const;
	int status(int i) const;

	const std::string get_contact_rate_map(int from, int to) const;
	double contact_rate(int from, int to) const;

	//===== AUGMENTED =====//

	//----- Augmented : Setter -----//
	void assign_augm(int i, const std::string &augmID = UNDEFINED, double value = std::numeric_limits<double>::quiet_NaN(),
					 const std::string &model_name = UNDEFINED, dictargs_ptr model_args_ptr = {}, dictargs range = {},
					 const std::string &rule_name = UNDEFINED, dictargs rule_args = {}, double rule_from = INF_DOUBLE);
	void set_augm_value(int i, const std::string &augmID, double x);	
	void augm_model_rand(int i, const std::string &augmID, sys_rnd_gen &prng);
	bool propose_valid_augm_move(int i, const std::string &augmID, sys_rnd_gen &prng);
	void augm_restore_value(int i, const std::string &augmID);
	void clear_augm(int i, const std::string &augmID);

	//----- Augmented : Getter -----//
	double augm(int i, const std::string &augmID) const;
	double augm_model_logpdf(int i, const std::string &augmID) const;
	double augm_model_logratio(int i, const std::string &augmID) const;
	double augm_propose_logratio(int i, const std::string &augmID) const;
	double augm_move_logratio(int i, const std::string &augmID) const;

	//===== LOGLIKELIHOOD =====//

	//----- Loglikelihood : Setter -----//
	void set_loglike(double x);
	void restore_loglike();

	//----- Loglikelihood : Getter -----//
	double get_loglike() const;
	double get_delta_loglike() const;

	//===== DIAGNOSTICS =====//
	void init_diagnostics(ulong sample_size);
	void update_diagnostics();
	std::vector<double> get_dgn_loglike() const;

private:
	//===== Pointer to Implementation =====//
	std::unique_ptr<Household_Impl> pImpl;
};

#endif // HOUSEHOLD__H