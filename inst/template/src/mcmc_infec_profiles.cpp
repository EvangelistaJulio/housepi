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

#include "mcmc_build.h"

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
			"Impossible to interpret value of \"" + val + "\" for \"" + profile->show_name() + "\". Try to provide a type.");

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

/* Store values of the infectivity profile. */
void MCMC::Impl::InfecProfile::init_lookup()
{
	if (profile->is_discrete())
	{
		bool is_suitable = true;
		// UPDATE : check if requires optional arguments
		try
		{
			UNUSED(__default_inst(profile->get_min() + 1));
		}
		catch (std::out_of_range &e)
		{
			is_suitable = false;
		}

		if (is_suitable)
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
		double x = __default_inst(i);
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
			tools::check_for_unused_entries(tools::headers_from_deflist(deflist),
											_headers);
			bool is_discrete = true;

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

			// if discrete and suitable initializes the lookup table
			init_lookup();
		}
	}
	catch (const std::exception &e)
	{
		tools::error({"While building infectivity profile definition:",
					  std::string(e.what())});
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  MCMC : Impl : InfecProfile : getter methods  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if a lookup table exists. */
bool MCMC::Impl::InfecProfile::is_lookup() const { return !_lookup_inst.empty(); }

/* Default method for instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::__default_inst(double t) const
{
	return profile->pdf(t);
}
/* Default method for cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::__default_cumul(double t) const
{
	return profile->cdf(t);
}

/* Lookup method for discrete instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::__lookup_inst(double t) const
{
	if ((t > 0) && (t < _size_lut))
		return _lookup_inst.at(static_cast<size_t>(t));
	else
		return 0;
}
/* Lookup method for discrete cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::__lookup_cumul(double t) const
{
	if (t < 0)
		return 0.;
	else if (t >= _size_lut)
		return 1.;
	else
		return _lookup_cumul.at(static_cast<size_t>(t));
}

/* Value of the instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::inst(double t) const
{
	if (profile->valid(t))
		return (this->*_do_inst)(t);
	else
		return 0.;
}
/* Value of the cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::cumul(double t) const
{
	if (profile->below(t))
		return 0.;
	else if (profile->above(t))
		return 1.;
	else
		return (this->*_do_cumul)(t);
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
void MCMC::Impl::update_infec_profile()
{
	if (_infec_profile.is_lookup())
		_infec_profile.update_lookup();
}

/* If infectivity profile is a lookup table, recover its buffered values. */
void MCMC::Impl::recover_infec_profile()
{
	if (_infec_profile.is_lookup())
		_infec_profile.recover_lookup();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : infectivity profile  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Virtual method for user initialization of infectivity profile. */
void MCMC::init_infec_profile() {};

/* Assign default infectivity profile during user model initialization. */
void MCMC::assign_template_infec_profile(std::map<std::string, std::string> x)
{
	pImpl->_infec_profile.assign_template(x);
}

/* Value of the instantaneous infectivity profile. */
double MCMC::infec_profile_inst(double t) const
{
	return pImpl->_infec_profile.inst(t);
}

/* Value of the cumulative infectivity profile. */
double MCMC::infec_profile_cumul(double t) const
{
	return pImpl->_infec_profile.cumul(t);
}
