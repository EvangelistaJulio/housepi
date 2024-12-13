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
	Household(const std::string &ID, size_t index);
	Household(Household &&other) noexcept;			  // Move constructor
	Household &operator=(Household &&other) noexcept; // Move assignment
	~Household();

	//===== SETTER =====//
	void assign_individual(const std::string &ID,
						   int status,
						   const dataIDmap<double> &data);
	bool change_status(indiv i, int new_status);
	void assign_contact_rate(indiv from, indiv to, const double *x);

	void disable_assert();

	//===== GETTER =====//
	const std::string &house_ID() const;
	size_t index() const;
	ushort size() const;
	ushort n_infectious() const;
	ushort n_susceptible() const;
	ushort n_missing_status() const;
	const std::set<indiv> &infectious_set() const;
	const std::set<indiv> &susceptible_set() const;
	const std::set<indiv> &missing_status_set() const;

	//===== INDIVIDUAL =====//

	//----- Individual : Setter -----//
	void assign_all_group_rates(indiv i, const rateIDmapvec<const double *> &x);

	//----- Individual : Getter -----//
	const std::string &ID(indiv i) const;
	int status(indiv i) const;
	double data(indiv i, const dataID &d) const;
	const double *data_pointer(indiv i, const dataID &d) const;
	double group_rate(indiv i, const rateID &r) const;

	double contact_rate(indiv from, indiv to) const;

	//===== AUGMENTED =====//

	//----- Augmented : Setter -----//
	void assign_augm(indiv i,
					 const augmID &a = UNDEFINED,
					 double value = std::numeric_limits<double>::quiet_NaN(),
					 bool is_discr = false,
					 const std::string &model_name = UNDEFINED,
					 dict_args_ptr model_args_ptr = {},
					 dict_args range = {},
					 const std::string &rule_name = UNDEFINED,
					 double rule_scale = INF_DOUBLE,
					 double rule_from = INF_DOUBLE);
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

	//===== EXTRA ATTRIBUTES =====//
	void _buffer_all_extra();
	void _restore_all_extra();

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
	const std::vector<double> &extra_vector(const extraID &e) const;

	//----- Matrix of doubles -----//
	void allocate_extra_matrix(const std::vector<extraID> &names, double xinit = 0.0);
	void set_extra_matrix(const extraID &e, size_t from, size_t to, double xnew);
	double extra_matrix(const extraID &e, size_t from, size_t to) const;
	const matrix<double> &extra_matrix(const extraID &e) const;

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
