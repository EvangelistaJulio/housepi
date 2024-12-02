/*
 * household.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__HOUSEHOLD__H
#define HOUSEPI__HOUSEHOLD__H

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

/* Forward declaration of Implementation. */
class Household_Impl;

class Household
{
public:
	//===== CONSTRUCTOR =====//
	Household();
	Household(const std::string &ID, int index);
	Household(Household &&other) noexcept;			  // Move constructor
	Household &operator=(Household &&other) noexcept; // Move assignment
	~Household();

	//===== SETTER =====//
	void assign_individual(const std::string &ID, int status, const dictargs &data);
	bool change_status(indiv i, int new_status);
	void init_contact_matrix();
	void assign_contact_rate(indiv from, indiv to, const std::string &parID);
	void set_contact_rate(indiv from, indiv to, const double *x);
	void disable_assert();

	//===== GETTER =====//
	const std::string &house_ID() const;
	int index() const;
	ushort size() const;
	ushort n_infectious() const;
	ushort n_susceptible() const;
	ushort n_missing_status() const;
	const std::set<indiv> &infectious_set() const;
	const std::set<indiv> &susceptible_set() const;
	const std::set<indiv> &missing_status_set() const;

	//===== INDIVIDUAL =====//

	//----- Individual : Setter -----//
	void assign_group_rate(indiv i, const std::string &rateID, const std::string &parID);
	void set_rate(indiv i, const std::string &rateID, const double *ptr);

	//----- Individual : Getter -----//
	const std::string &ID(indiv i) const;
	double data(indiv i, const std::string &dataID) const;
	const double *data_pointer(indiv i, const std::string &dataID) const;
	const std::unordered_map<std::string, std::string> &get_group_rate_map(indiv i) const;
	double group_rate(indiv i, const std::string &rateID) const;
	int status(indiv i) const;

	const std::string &get_contact_rate_map(indiv from, indiv to) const;
	double contact_rate(indiv from, indiv to) const;

	//===== AUGMENTED =====//

	//----- Augmented : Setter -----//
	void assign_augm(indiv i, const std::string &augmID = UNDEFINED, double value = std::numeric_limits<double>::quiet_NaN(),
					 const std::string &model_name = UNDEFINED, dictargs_ptr model_args_ptr = {}, dictargs range = {},
					 const std::string &rule_name = UNDEFINED, dictargs rule_args = {}, double rule_from = INF_DOUBLE);
	void set_augm_value(indiv i, const std::string &augmID, double x);
	void augm_model_rand(indiv i, const std::string &augmID, sys_rnd_gen &prng);
	bool propose_valid_augm_move(indiv i, const std::string &augmID, sys_rnd_gen &prng);
	void restore_augm_value(indiv i, const std::string &augmID);
	void clear_augm(indiv i, const std::string &augmID);

	//----- Augmented : Getter -----//
	double augm(indiv i, const std::string &augmID) const;
	double augm_model_logpdf(indiv i, const std::string &augmID) const;
	double augm_model_logratio(indiv i, const std::string &augmID) const;
	double augm_propose_logratio(indiv i, const std::string &augmID) const;
	double augm_move_logratio(indiv i, const std::string &augmID) const;

	//===== EXTRA ATTRIBUTES =====//
	void _buffer_all_extra();
	void _restore_all_extra();

	void allocate_extra_int(const std::vector<std::string> names, int init_value = 0.0);
	void set_extra_int(const std::string &intID, int new_value);
	int extra_int(const std::string &intID) const;

	void allocate_extra_double(const std::vector<std::string> names, double init_value = 0.0);
	void set_extra_double(const std::string &dblID, double new_value);
	double extra_double(const std::string &dblID) const;

	void allocate_extra_vector(const std::vector<std::string> names, double init_value = 0.0);
	void set_extra_vector(const std::string &vecID, size_t i, double new_value);
	double extra_vector(const std::string &vecID, size_t i) const;
	const std::vector<double> &extra_vector(const std::string &vecID) const;

	void allocate_extra_matrix(const std::vector<std::string> names, double init_value = 0.0);
	void set_extra_matrix(const std::string &mtxID, size_t from, size_t to, double new_value);
	double extra_matrix(const std::string &mtxID, size_t from, size_t to) const;
	const squarematrix &extra_matrix(const std::string &mtxID) const;

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
	const std::vector<double> &get_dgn_loglike() const;

private:
	//===== Pointer to Implementation =====//
	std::unique_ptr<Household_Impl> pImpl;
};

#endif // HOUSEPI__HOUSEHOLD__H
