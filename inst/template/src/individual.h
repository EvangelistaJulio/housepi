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
/////////////////////////////////////  Group Rate : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class GroupRate
{
public:
	std::string name = "";

	GroupRate();
	GroupRate(const std::string &name,
			  const std::vector<const double *> &pointers);
	~GroupRate();

	double get() const;

private:
	std::vector<const double *> pointers = {};

	double (GroupRate::*_func_getter_)() const = &GroupRate::_getter_error_;
	double _getter_error_() const;
	double _getter_default_() const;
	double _getter_multiple_() const;
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

/* Class for Individual household study member. */
class Individual
{
public:
	//===== CONSTRUCTOR =====//
	Individual();
	Individual(const std::string &ID,
			   int status,
			   const dataIDmap<double> &data = {});
	~Individual();

	//===== SETTER =====//
	void set_status(int s);
	void disable_assert();
	void assign_all_group_rates(const rateIDmapvec<const double *> &group_rate_ptrs);

	//===== GETTER =====//
	const std::string &get_ID() const;
	int get_status() const;
	double get_data(const dataID &d) const;
	const double *get_data_pointer(const dataID &d) const;
	double get_group_rate(const rateID &r) const;

	//===== AUGMENTED VARIABLES =====//

	//---- Getter -----//
	double get_augm(const augmID &a) const;
	double augm_model_logpdf(const augmID &a) const;
	double augm_model_logratio(const augmID &a) const;
	double augm_propose_logratio(const augmID &a) const;
	double augm_move_logratio(const augmID &a) const;

	//---- Setter -----//
	void clear_augm(const augmID &a);
	void assign_augm(const augmID &a = UNDEFINED,
					 const std::string &name = UNDEFINED,
					 double value = std::numeric_limits<double>::quiet_NaN(),
					 bool is_discrete = false,
					 const std::string &model_name = UNDEFINED,
					 const dict_args_ptr &model_args_ptr = {},
					 const dict_args &range = {},
					 const std::string &rule_name = UNDEFINED,
					 double rule_scale = INF_DOUBLE,
					 double rule_from = INF_DOUBLE);
	void set_augm_value(const augmID &a, double x);
	void augm_model_rand(const augmID &a, rnd_gen &prng);
	bool propose_valid_augm_move(const augmID &a, rnd_gen &prng);
	void restore_augm_value(const augmID &a);

private:
	std::string ID = UNDEFINED;	 // user ID
	int status = 0;				 // status of the individual
	dataIDmap<double> data = {}; // database associated to individual

	rateIDmap<GroupRate> rate_ptr = {};
	augmIDmap<Base_Variable> augm_map = {};

	//===== Auxxiliary Assert =====//
	bool assert_data(const dataID &d) const;
	bool assert_augm(const augmID &a) const;

	//===== Setter (auxiliary methods) =====//
	void (Individual::*_func_set_augm_)(const augmID &, double) = &Individual::_set_augm_assert_;
	void _set_augm_assert_(const augmID &a, double xnew);
	void _set_augm_default_(const augmID &a, double xnew);

	//===== Getter (auxiliary methods) =====//
	double (Individual::*_func_get_augm_)(const augmID &) const = &Individual::_get_augm_assert_;
	double _get_augm_assert_(const augmID &a) const;
	double _get_augm_default_(const augmID &a) const;

	double (Individual::*_func_get_data_)(const dataID &) const = &Individual::_get_data_assert_;
	double _get_data_assert_(const dataID &d) const;
	double _get_data_default_(const dataID &d) const;

	const double *(Individual::*_func_get_data_pointer_)(const dataID &) const = &Individual::_get_data_pointer_assert_;
	const double *_get_data_pointer_assert_(const dataID &d) const;
	const double *_get_data_pointer_default_(const dataID &d) const;
};

#endif // HOUSEPI__INDIVIDUAL__H
