/*
 * individual.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__INDIVIDUAL__H
#define HOUSEPI__INDIVIDUAL__H

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

#include "variable.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Individual : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Class for Individual household study member. */
class Individual
{
public:
	//===== CONSTRUCTOR =====//
	Individual();
	Individual(const std::string &new_ID, int new_status, const dictargs &data);
	~Individual();

	//===== SETTER =====//
	void set_status(int s);
	void disable_assert();
	void assign_all_group_rates(const dict_rates_ptr &group_rate_ptrs);

	//===== GETTER =====//
	const std::string &get_ID() const;
	int get_status() const;
	double get_data(const dataID &d) const;
	const double *get_data_pointer(const dataID &d) const;
	double get_group_rate(const std::string &rateID) const;

	//===== AUGMENTED VARIABLES =====//

	//---- Getter -----//
	double get_augm(const augmID &a) const;
	double augm_model_logpdf(const augmID &a) const;
	double augm_model_logratio(const augmID &a) const;
	double augm_propose_logratio(const augmID &a) const;
	double augm_move_logratio(const augmID &a) const;

	//---- Setter -----//
	void assign_augm(const augmID &a = UNDEFINED,
					 double value = std::numeric_limits<double>::quiet_NaN(),
					 bool is_discr = false,
					 const std::string &model_name = UNDEFINED, 
					 dictargs_ptr model_args_ptr = {}, 
					 dictargs range = {},
					 const std::string &rule_name = UNDEFINED, 
					 double rule_scale = INF_DOUBLE, 
					 double rule_from = INF_DOUBLE);
	void set_augm_value(const augmID &a, double x);
	void augm_model_rand(const augmID &a, sys_rnd_gen &prng);
	bool propose_valid_augm_move(const augmID &a, sys_rnd_gen &prng);
	void restore_augm_value(const augmID &a);
	void clear_augm(const augmID &a);

private:
	std::string ID = UNDEFINED; // user ID
	int status = 0;				// status of the individual
	dictargs data = {};			// database associated to individual

	std::unordered_map<std::string, std::vector<const double *>> rate_ptr = {};
	std::unordered_map<std::string, Base_Variable> augm_map = {};

	//===== Auxxiliary Assert =====//
	bool assert_data(const dataID &d) const;
	bool assert_augm(const augmID &a) const;

	//===== Setter (auxiliary methods) =====//
	void (Individual::*_func_set_augm_)(const std::string &, double) = &Individual::_set_augm_assert_;
	void _set_augm_assert_(const augmID &a, double new_value);
	void _set_augm_default_(const augmID &a, double new_value);

	//===== Getter (auxiliary methods) =====//
	double (Individual::*_func_get_augm_)(const std::string &) const = &Individual::_get_augm_assert_;
	double _get_augm_assert_(const augmID &a) const;
	double _get_augm_default_(const augmID &a) const;

	double (Individual::*_func_get_data_)(const std::string &) const = &Individual::_get_data_assert_;
	double _get_data_assert_(const dataID &d) const;
	double _get_data_default_(const dataID &d) const;

	const double *(Individual::*_func_get_data_pointer_)(const std::string &) const = &Individual::_get_data_pointer_assert_;
	const double *_get_data_pointer_assert_(const dataID &d) const;
	const double *_get_data_pointer_default_(const dataID &d) const;
};

#endif // HOUSEPI__INDIVIDUAL__H
