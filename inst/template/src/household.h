/*
 * household.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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

class Household
{
public:
	//=======================//
	//===== CONSTRUCTOR =====//
	//=======================//

	Household(const std::string &ID = NULL_STRING,
			  size_t index = SIZE_MAX);
	Household(Household &&other) noexcept;			  // Move constructor
	Household &operator=(Household &&other) noexcept; // Move assignment
	~Household();

	//==================//
	//===== SETTER =====//
	//==================//

	void assign_individual(const std::string &ID,
						   int status,
						   const dataIDmap<double> &data);
	bool change_status(indiv i, int new_status);
	void assign_contact_rate(indiv from, indiv to, const double *x);

	//==================//
	//===== GETTER =====//
	//==================//

	const std::string &house_ID() const;
	size_t index() const;
	ushort size() const;
	ushort n_infectious() const;
	ushort n_susceptible() const;
	ushort n_missing_status() const;
	const std::set<indiv> &infectious_set() const;
	const std::set<indiv> &susceptible_set() const;
	const std::set<indiv> &missing_status_set() const;

	//======================//
	//===== INDIVIDUAL =====//
	//======================//

	//----- Individual : Setter -----//
	void assign_all_group_rates(indiv i, const rateIDmapvec<const double *> &x);

	//----- Individual : Getter -----//
	const std::string &ID(indiv i) const;
	int status(indiv i) const;
	double data(indiv i, const dataID &d) const;
	double group_rate(indiv i, const rateID &r) const;

	double contact_rate(indiv from, indiv to) const;

	//=====================//
	//===== AUGMENTED =====//
	//=====================//

	//----- Augmented : Setter -----//
	void assign_augm(indiv i,
					 const augmID &a = NULL_STRING,
					 const std::string &name = NULL_STRING,
					 double value = std::numeric_limits<double>::quiet_NaN(),
					 bool is_discrete = false,
					 const std::string &model_name = NULL_STRING,
					 const std::string &rule_name = NULL_STRING,
					 double rule_scale = NULL_DOUBLE,
					 double rule_from = NULL_DOUBLE);

	void assign_augm_model_scalar_arg(indiv i, const augmID &a,
									  const std::string &arg_name,
									  double arg_value);
	void assign_augm_model_pointer_arg(indiv i, const augmID &a,
									   const std::string &arg_name,
									   const double *ptr_to_arg_value);

	void assign_augm_model_scalar_range_min(indiv i, const augmID &a,
											double arg_value);
	void assign_augm_model_pointer_range_min(indiv i, const augmID &a,
											 const double *ptr_to_arg_value);

	void assign_augm_model_scalar_range_max(indiv i, const augmID &a,
											double arg_value);
	void assign_augm_model_pointer_range_max(indiv i, const augmID &a,
											 const double *ptr_to_arg_value);
	void init_augm_model_reshaping(indiv i, const augmID &a);

	void set_augm_value(indiv i, const augmID &a, double x);
	void augm_model_rand(indiv i, const augmID &a, rnd_gen &prng);
	bool propose_valid_augm_move(indiv i, const augmID &a, rnd_gen &prng);
	void restore_augm_value(indiv i, const augmID &a);
	void clear_augm(indiv i, const augmID &a);

	//----- Augmented : Getter -----//
	double augm(indiv i, const augmID &a) const;
	double augm_model_logpdf(indiv i, const augmID &a) const;
	double augm_model_logratio(indiv i, const augmID &a) const;
	double augm_propose_logratio(indiv i, const augmID &a) const;
	double augm_move_logratio(indiv i, const augmID &a) const;

	//============================//
	//===== EXTRA ATTRIBUTES =====//
	//============================//

	// FIXME: am I buffering right if infectivity profile depends on param to infer?
	void buffer_all_extra();

	//----- Single integers -----//
	void allocate_extra_integer(const std::vector<extraID> &names, int xinit = 0.0);
	void set_extra_integer(const extraID &e, int xnew);
	int extra_integer(const extraID &e) const;

	//----- Single doubles -----//
	void allocate_extra_double(const std::vector<extraID> &names, double xinit = 0.0);
	void set_extra_double(const extraID &e, double xnew);
	double extra_double(const extraID &e) const;

	//----- Vector of doubles -----//
	void allocate_extra_vector(const std::vector<extraID> &names, double xinit = 0.0);
	void set_extra_vector(const extraID &e, size_t i, double xnew);
	double extra_vector(const extraID &e, size_t i) const;
	const dbl_vec &extra_vector(const extraID &e) const;

	//----- Matrix of doubles -----//
	void allocate_extra_matrix(const std::vector<extraID> &names, double xinit = 0.0);
	void set_extra_matrix(const extraID &e, size_t from, size_t to, double xnew);
	double extra_matrix(const extraID &e, size_t from, size_t to) const;
	const matrix<double> &extra_matrix(const extraID &e) const;

	//=========================//
	//===== LOGLIKELIHOOD =====//
	//=========================//

	//----- Loglikelihood : Setter -----//
	void set_loglike(double x);
	void restore_loglike();

	//----- Loglikelihood : Getter -----//
	double get_loglike() const;
	double get_delta_loglike() const;

	//=======================//
	//===== DIAGNOSTICS =====//
	//=======================//

	void init_diagnostics(ulong sample_size);
	void update_diagnostics(ulong count);
	const dbl_vec &get_dgn_loglike() const;

private:
	//===== Pointer to Implementation =====//
	class Impl;
	std::unique_ptr<Impl> pImpl;
};

#endif // HOUSEPI__HOUSEHOLD__H
