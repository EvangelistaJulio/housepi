/*
 * profiles.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "distribution.h"

// UPDATE: this should be merge into functions finally.
// It could be done by changing Arguments.

/* Base Class for Infectivity Profiles.
	NOTE : this classes access to parameters passing through MCMC class.
	This way I can check which parameters changes augmented dependencies. */
class Base_InfecProfile
{
public:
	//===== CONSTRUCTOR =====//
	Base_InfecProfile(std::function<double(parID)> get = nullptr);
	virtual ~Base_InfecProfile();

	//===== GETTER =====//
	double get_arg(ushort i);

	virtual double pdf(const std::vector<double> &x);
	virtual double cdf(const std::vector<double> &x);

protected:
	//===== SETTER =====//
	void _init_defaults_(const std::string &name,
						 const std::vector<std::string> &arg_name,
						 const std::vector<double> &arg_value);

private:
	std::string name = UNDEFINED;
	str_vec _arg_names_ = {};
	std::vector<double> _values_ = {};
	std::vector<parID> _indexes_ = {};

	//----- Functional methods -----//
	double (Base_InfecProfile::*_func_get_arg_)(ushort) = &Base_InfecProfile::_internal_get_arg_;
	std::function<double(parID)> _external_param_getter_ = nullptr;

	double _internal_get_arg_(ushort i);
	double _external_get_arg_(ushort i);
};

/* Assign the infectivity profile. */
std::unique_ptr<Base_InfecProfile> infecprofile_switchboard(const std::string &name,
															std::function<double(parID)> get = nullptr);
