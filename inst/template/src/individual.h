/*
 * individual.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef INDIVIDUAL__H
#define INDIVIDUAL__H

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

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
	void set_rate(const std::string &rateID, const double *x);	
	void assign_group_rate(const std::string &rateID, const std::string &parID);

	//===== GETTER =====//
	std::string get_ID() const;
	int get_status() const;
	double get_data(const std::string &dataID) const;
	const double *get_data_pointer(const std::string &dataID) const;
	double get_group_rate(const std::string &rateID) const;
	std::unordered_map<std::string, std::string> get_group_rate_map() const;

	//===== AUGMENTED VARIABLES =====//

	//---- Getter -----//
	double get_augm(const std::string &augmID) const;
	double augm_model_logpdf(const std::string &augmID) const;
	double augm_model_logratio(const std::string &augmID) const;
	double augm_propose_logratio(const std::string &augmID) const;
	double augm_move_logratio(const std::string &augmID) const;

	//---- Setter -----//
	void assign_augm(const std::string &augmID = UNDEFINED, double value = std::numeric_limits<double>::quiet_NaN(),
					 const std::string &model_name = UNDEFINED, dictargs_ptr model_args_ptr = {}, dictargs range = {},
					 const std::string &rule_name = UNDEFINED, dictargs rule_args = {}, double rule_from = INF_DOUBLE);
	void set_augm_value(const std::string &augmID, double x);
	void augm_model_rand(const std::string &augmID, sys_rnd_gen &prng);
	bool propose_valid_augm_move(const std::string &augmID, sys_rnd_gen &prng);
	void augm_restore_value(const std::string &augmID);
	void clear_augm(const std::string &augmID);

private:
	std::string ID = UNDEFINED; // user ID
	int status = 0;				// status of the individual
	dictargs data = {};			// database associated to individual

	dictargs_ptr rate_ptr = {};
	std::unordered_map<std::string, std::string> _which_rate_ = {};
	std::unordered_map<std::string, Base_Variable> augm_var = {};

	//===== Auxxiliary Assert =====//
	void assert_data(const std::string &dataID) const;
	void assert_augm(const std::string &augmID) const;

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

#endif // INDIVIDUAL__H