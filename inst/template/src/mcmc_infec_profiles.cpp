/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_infec_profiles.cpp
 *
 */

#include "mcmc_impl.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC : Impl : assign infectivity profile arguments  ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Adjust lists of parameters to be inferred making parameter as dependent on cache. */
void MCMC::Impl::make_param_cache_dependent(const std::string &par_name)
{
	// remove a parameter from standard if it was defined as such
	auto it = std::find(param_standard_infer_list.begin(), param_standard_infer_list.end(), par_name);
	if (it != param_standard_infer_list.end())
		param_standard_infer_list.erase(it);

	// add the parameter as cache
	it = std::find(param_cache_infer_list.begin(), param_cache_infer_list.end(), par_name);
	if (it == param_cache_infer_list.end())
		param_cache_infer_list.push_back(par_name);
}

/* Assign arguments to infectivity profile. */
void MCMC::Impl::InfecProfile::_assign_arguments(const str_vec &keys,
												 const str_vec &vals,
												 const str_vec &typs)
{
	// check size validity of keys, values and type
	if ((keys.size() != vals.size()) || (typs.size() != vals.size()))
		throw std::invalid_argument(
			"Mismatching definitions in arguments.");

	// loop over all arguments
	for (size_t k = 0; k < keys.size(); ++k)
		_assign_single_attribute(keys[k], vals[k], typs[k], ATTR_ARGUMENT);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : InfecProfile : assign range  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign general attribute of infectivity profile. */
void MCMC::Impl::InfecProfile::_assign_single_attribute(const std::string &key,
														const std::string &val,
														const std::string &at,
														bool is_arg)
{
	if (tools::is_na(val) && tools::is_na(at))
		; // pass

	else if (!tools::is_na(val) && tools::is_na(at))
		// UPDATE: this error is not user friendly sinne API is supposed to associated type
		throw std::invalid_argument(
			"Unable to interpret value of \"" + val + "\" for \"" + profile->show_name() + "\". Try to provide a type.");

	else if (tools::is_na(val) && !tools::is_na(at))
		throw std::invalid_argument(
			"Found type without a paired value for \"" + profile->show_name() + "\".");

	else if (at == TYPE_PARAM)
	{
		if (!pImpl.param_map.at(val).is_constant())
		{
			profile->assign_attr(key, is_arg, pImpl.get_param_pointer(val));
			pImpl.make_param_cache_dependent(val);
		}
		else
			profile->assign_attr(key, is_arg, pImpl.get_param(val));
	}

	else if (at == TYPE_SCALAR)
		profile->assign_attr(key, is_arg, tools::stod(val, "argument with scalar type"));

	else
		throw std::invalid_argument(
			"Unsopported type \"" + at + "\" found while loading \"" + val + "\" for \"" + profile->show_name() + "\".");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : Impl : InfecProfile : assign template  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Function to assign a template to infectivity profile. */
void MCMC::Impl::InfecProfile::assign_template(std::map<std::string, std::string> &x)
{
	if (deflist.empty())
	{
		if (!x.count("name"))
			throw std::invalid_argument("Missing name for infectivity profile template.");

		for (const std::string &key : _headers)
			if (!x.count(key))
				x.emplace(key, HOUSEPI_NA_STRING);

		if (tools::is_na(x.at("discr")))
			x["discr"] = "false";

		deflist = x;
		build();
	}
}

/* Append infectivity profile to augmented variable definitions. */
void MCMC::Impl::InfecProfile::append_augm_definition()
{
	if (is_augmented())
		for (const auto &pair : _augm_infec_profile_deflist)
		{
			auto &aux = pImpl._augm_deftable.at(pair.first);
			aux.insert(aux.end(), pair.second.begin(), pair.second.end());
		}
}

/* Store values of the infectivity profile. */
void MCMC::Impl::InfecProfile::init_lookup()
{
	if (profile->is_discrete())
	{
		if (tools::is_na(profile->get_max()))
			_size_lut = DEFAULT_MAX_SIZE;
		else
			_size_lut = profile->get_max() + 1;

		_lookup_inst.resize(_size_lut);
		_lookup_cumul.resize(_size_lut);
		_lookup_buffer_inst.resize(_size_lut);
		_lookup_buffer_cumul.resize(_size_lut);

		update_lookup();

		_do_inst = &InfecProfile::__lookup_inst;
		_do_cumul = &InfecProfile::__lookup_cumul;
	}
	else
	{
		_do_inst = &InfecProfile::__default_inst;
		_do_cumul = &InfecProfile::__default_cumul;
	}
}

/* Restore buffered lookup table. */
void MCMC::Impl::InfecProfile::recover_lookup()
{
	_lookup_inst.swap(_lookup_buffer_inst);
	_lookup_cumul.swap(_lookup_buffer_cumul);
}

/* Store values of the infectivity profile. */
void MCMC::Impl::InfecProfile::update_lookup()
{
	recover_lookup();

	// Set infectivity to 0 on the same day
	_lookup_inst[0] = 0.;
	_lookup_cumul[0] = 0.;

	for (size_t i = 1; i < _size_lut; ++i)
	{
		Individual dummy;
		double x = __default_inst(dummy, i);
		_lookup_inst[i] = x;
		_lookup_cumul[i] = x + _lookup_cumul[i - 1];
	}
}

/* Build infectivity profile distribution . */
void MCMC::Impl::InfecProfile::build()
{
	try
	{
		if (!deflist.empty())
		{
			// check entries
			tools::check_for_unused_entries(
				tools::headers_from_deflist(deflist), _headers);

			// check if definition is generic or augmented
			str_vec types = tools::parse_list_args(deflist.at("arg_at"));
			types += {deflist.at("min_at"), deflist.at("max_at")};

			bool is_augm = false;
			for (const std::string &at : types)
				if (at == TYPE_AUGM)
					is_augm = true;

			if (is_augm)
				_build_augmented();
			else
				_build_generic();
		}
	}
	catch (const std::exception &e)
	{
		tools::error("While building infectivity profile:\n" +
					 std::string(e.what()));
	}
}

/* Build for generic infectivity profile. */
void MCMC::Impl::InfecProfile::_build_generic()
{
	bool is_discrete = true;
	_augm_infec_profile_deflist.clear();

	//----- Name -----//
	if ((deflist.count("name")) && (!tools::is_na(deflist.at("name"))))
	{
		profile = CompositeDistrRegistry(
					  "infectivity profile",
					  {std::make_shared<InfecProfileRegistry>(),
					   std::make_shared<ContinuosDistrRegistry>(),
					   std::make_shared<DiscreteDistrRegistry>()})
					  .get_distr(deflist.at("name"));
	}
	else
		throw std::invalid_argument("Missing name.");

	//----- Discrete -----//
	if (!profile->is_discrete()) // not defined discrete
		if (deflist.count("discr"))
			is_discrete = tools::stobool(deflist.at("discr"));

	//----- Arguments -----//
	if (deflist.count("arg_key") && deflist.count("arg_val") && deflist.count("arg_at"))
		_assign_arguments(tools::parse_list_args(deflist.at("arg_key")),
						  tools::parse_list_args(deflist.at("arg_val")),
						  tools::parse_list_args(deflist.at("arg_at")));

	//----- Range minimum -----//
	_assign_single_attribute("min",
							 (is_discrete) ? "1" : "0", TYPE_SCALAR,
							 ATTR_RANGE);

	if (deflist.count("min") && deflist.count("min_at"))
		_assign_single_attribute("min",
								 deflist.at("min"), deflist.at("min_at"),
								 ATTR_RANGE);

	// check infectivity profile minimum is correct
	if (is_discrete)
	{
		if (profile->get_min() < 1.)
			throw std::invalid_argument(
				"Minimum of discrete infectivity profile must be >= 1.");
	}
	else if (profile->get_min() < 0.)
		throw std::invalid_argument(
			"Minimum of continuos infectivity profile must be >= 0.");

	//----- Range maximum -----//
	if (deflist.count("max") && deflist.count("max_at"))
		_assign_single_attribute("max",
								 deflist.at("max"), deflist.at("max_at"),
								 ATTR_RANGE);

	//----- Reshaping -----//
	profile->init_reshaping(is_discrete);

	// if discrete and generic, initializes lookup table
	init_lookup();
}

/* Build for individual augmented function infectivity profile. */
void MCMC::Impl::InfecProfile::_build_augmented()
{
	std::string func_name;

	profile = nullptr;
	_do_inst = &InfecProfile::__error;
	_do_cumul = &InfecProfile::__error;

	if ((deflist.count("name")) && (!tools::is_na(deflist.at("name"))))
		func_name = deflist.at("name");
	else
		throw std::invalid_argument("Missing name.");

	// initialize for generation time, instantaneous and cumulative profiles
	for (auto &entry : pImpl.augm_file_headers)
		_augm_infec_profile_deflist.emplace(entry, str_vec(3, HOUSEPI_NA_STRING));

	_augm_infec_profile_deflist.at("const") = str_vec({tools::truth(true)}) + str_vec(2, tools::truth(false));
	_augm_infec_profile_deflist.at("init_rand") = str_vec(3, tools::truth(false));

	// match infectivity profile and augmented variable definition
	const dict_str_vec matcher = {
		{"name", str_vec({AUX_IP_GEN_TIME, AUX_IP_INST_PROFILE, AUX_IP_CUMUL_PROFILE})},
		{"func", str_vec({HOUSEPI_NA_STRING, "d" + func_name, "p" + func_name})},
		{"farg_key", str_vec({HOUSEPI_NA_STRING}) + str_vec(2, deflist.at("arg_key"))},
		{"farg_val", str_vec({HOUSEPI_NA_STRING}) + str_vec(2, deflist.at("arg_val"))},
		{"farg_at", str_vec({HOUSEPI_NA_STRING}) + str_vec(2, deflist.at("arg_at"))},
		{"discr", str_vec(3, deflist.at("discr"))},
		{"min", str_vec({HOUSEPI_NA_STRING}) + str_vec(2, deflist.at("min"))},
		{"min_at", str_vec({HOUSEPI_NA_STRING}) + str_vec(2, deflist.at("min_at"))},
		{"max", str_vec({HOUSEPI_NA_STRING}) + str_vec(2, deflist.at("max"))},
		{"max_at", str_vec({HOUSEPI_NA_STRING}) + str_vec(2, deflist.at("max_at"))}};

	for (auto &pair : matcher)
		_augm_infec_profile_deflist.at(pair.first) = pair.second;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  MCMC : Impl : InfecProfile : getter methods  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Default method for instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::__default_inst(Individual &, double t)
{
	if (profile->valid(t))
		return profile->pdf(t);
	else
		return 0.;
}
/* Default method for cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::__default_cumul(Individual &, double t)
{
	if (profile->below(t))
		return 0.;
	else if (profile->above(t))
		return 1.;
	else
		return profile->cdf(t);
}

/* Lookup method for discrete instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::__lookup_inst(Individual &, double t)
{
	if (profile->valid(t) && (t < _size_lut))
		return _lookup_inst.at(static_cast<size_t>(t));
	else
		return 0;
}
/* Lookup method for discrete cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::__lookup_cumul(Individual &, double t)
{
	if (profile->below(t))
		return 0.;
	else if ((t >= _size_lut) || (profile->above(t)))
		return 1.;
	else
		return _lookup_cumul.at(static_cast<size_t>(t));
}

/* Value of the instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::__augm_inst(Individual &member, double t)
{
	member.get_augm_variable(AUX_IP_GEN_TIME).set_value(t);
	return member.augm(AUX_IP_INST_PROFILE);
}

/* Value of the cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::__augm_cumul(Individual &member, double t)
{
	member.get_augm_variable(AUX_IP_GEN_TIME).set_value(t);
	return member.augm(AUX_IP_CUMUL_PROFILE);
}

/* Error method. */
double MCMC::Impl::InfecProfile::__error(Individual &, double)
{
	throw std::runtime_error(
		"Caled generic infectivity profile while augmented for each individual.");
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : update\recover lookup tables  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* If infectivity profile is a lookup table, store its values. */
void MCMC::Impl::_update_infec_profile()
{
	if (_infec_profile.is_lookup())
		_infec_profile.update_lookup();
}

/* If infectivity profile is a lookup table, recover its buffered values. */
void MCMC::Impl::recover_infec_profile()
{
	// UPDATE: this should not be recoverd if constant for example...
	if (_infec_profile.is_lookup())
		_infec_profile.recover_lookup();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : infectivity profile  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign default infectivity profile during user model initialization. */
void MCMC::assign_template_infec_profile(std::map<std::string, std::string> x)
{
	pImpl->_infec_profile.assign_template(x);
}

/* Value of the instantaneous infectivity profile. */
double MCMC::infec_profile_inst(Individual &infector, double t) const
{
	return pImpl->_infec_profile.inst(infector, t);
}

/* Value of the cumulative infectivity profile. */
double MCMC::infec_profile_cumul(Individual &infector, double t) const
{
	return pImpl->_infec_profile.cumul(infector, t);
}
