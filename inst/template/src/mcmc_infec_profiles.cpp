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

/* Assign arguments to infectivity profile. */
void MCMC::Impl::ImplInfecProfile::assign_arguments_(const str_vec &keys,
													 const str_vec &vals,
													 const str_vec &typs)
{
	// check size validity of keys, values and type
	if ((keys.size() != vals.size()) || (typs.size() != vals.size()))
		throw std::invalid_argument(
			"Mismatching definitions in arguments.");

	// loop over all arguments
	for (size_t k = 0; k < keys.size(); ++k)
		assign_single_attribute_(keys[k], vals[k], typs[k], HSP_ATTR_ARGUMENT);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : InfecProfile : assign range  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign general attribute of infectivity profile. */
void MCMC::Impl::ImplInfecProfile::assign_single_attribute_(const std::string &key,
															const std::string &val,
															const std::string &at,
															bool is_arg)
{
	if (tools::is_na(val) && tools::is_na(at))
		; // pass

	else if (!tools::is_na(val) && tools::is_na(at))
		// UPDATE: this error is not user friendly sinne API is supposed to associated type
		throw std::invalid_argument(
			"Unable to interpret value of \"" + val + "\" for \"" + profile_->show_name() + "\". Try to provide a type.");

	else if (tools::is_na(val) && !tools::is_na(at))
		throw std::invalid_argument(
			"Found type without a paired value for \"" + profile_->show_name() + "\".");

	else if (tools::stoupper(at) == TYPE_PARAM)
	{
		if (pImpl_.param.dictionary.at(val).is_constant())
			profile_->assign_attr(key, is_arg, pImpl_.param.get(val));
		else
		{
			profile_->assign_attr(key, is_arg, pImpl_.param.get_ptr(val));
			pImpl_.param.make_cache_dependent(val);
		}
	}

	else if (tools::stoupper(at) == TYPE_SCALAR)
		profile_->assign_attr(key, is_arg, tools::stod(val, "argument with scalar type"));

	else
		throw std::invalid_argument(
			"Unsopported type \"" + at + "\" found while loading \"" + val + "\" for \"" + profile_->show_name() + "\".");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : Impl : InfecProfile : assign template  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign model template. */
void MCMC::Impl::ImplInfecProfile::assign_template(const dict_str_str &deflist) { deflist_template_ = deflist; }

/* Initialize infectivity profile with template if custom not provied. */
void MCMC::Impl::ImplInfecProfile::init()
{
	// if not previous build proceed with template
	if ((profile_ == nullptr) && augm_infec_profile_deflist_.empty())
		// if no template, assume no need for an infectivity profile
		if (!deflist_template_.empty())
		{
			if (!deflist_template_.count("name"))
				throw std::invalid_argument(
					"Missing name for infectivity profile template.");

			for (const std::string &key : headers_)
				if (!deflist_template_.count(key))
					deflist_template_.emplace(key, HSP_NA_STRING);

			std::string str_is_discr = "false";
			if (deflist_template_.count("discr"))
				if (!tools::is_na(deflist_template_.at("discr")))
					str_is_discr = deflist_template_.at("discr");
			deflist_template_["discr"] = str_is_discr;

			build(deflist_template_);
		}
}

/* Append infectivity profile to augmented variable definitions.
 * UPDATE: problematic in case of redefinitions of infectivity profile
 */
void MCMC::Impl::ImplInfecProfile::append_augm_definition()
{
	if (is_augmented())
	{
		// check for reserved augmented names
		for (const std::string &a : pImpl_.augm.deftable.at("name"))
			for (const std::string &forbid : {HSP_PROFILE_INST,
											  HSP_PROFILE_CUMUL})
				if (a == forbid)
					throw std::invalid_argument(
						"Augmented variable \"" + forbid + "\" is reserved. Please, pick another name.");

		// append new definitions
		for (const auto &pair : augm_infec_profile_deflist_)
		{
			auto &aux = pImpl_.augm.deftable.at(pair.first);
			aux.insert(aux.end(), pair.second.begin(), pair.second.end());
		}
	}
}

/* Store values of the infectivity profile. */
void MCMC::Impl::ImplInfecProfile::init_lookup_()
{
	if (profile_->is_discrete())
	{
		inst_ = &ImplInfecProfile::do_inst_lookup_;
		cumul_ = &ImplInfecProfile::do_cumul_lookup_;

		if (tools::is_na(profile_->get_max()))
			size_lut_ = DEFAULT_MAX_SIZE;
		else
			size_lut_ = profile_->get_max() + 1;

		lookup_inst_.resize(size_lut_);
		lookup_cumul_.resize(size_lut_);
		lookup_buffer_inst_.resize(size_lut_);
		lookup_buffer_cumul_.resize(size_lut_);

		update(true);
	}
}

/* If infectivity profile is a lookup table, store its values. */
void MCMC::Impl::ImplInfecProfile::update(bool is_cache_depen)
{
	// UPDATE: this should not be updated if constant
	if (is_cache_depen)
		if (is_lookup())
		{
			// move current lookup table to buffer
			lookup_inst_.swap(lookup_buffer_inst_);
			lookup_cumul_.swap(lookup_buffer_cumul_);

			// Set infectivity to 0 on the same day
			lookup_inst_[0] = 0.;
			lookup_cumul_[0] = 0.;

			for (size_t i = 1; i < size_lut_; ++i)
			{
				double x = do_inst_default_(nullptr, i);
				lookup_inst_[i] = x;
				lookup_cumul_[i] = x + lookup_cumul_[i - 1];
			}
		}
}

/* If infectivity profile is a lookup table, recover its buffered values. */
void MCMC::Impl::ImplInfecProfile::recover(bool is_cache_depen)
{
	// UPDATE: this should not be recoverd if constant
	if (is_cache_depen)
		if (is_lookup())
		{
			lookup_inst_.swap(lookup_buffer_inst_);
			lookup_cumul_.swap(lookup_buffer_cumul_);
		}
}

/* Build infectivity profile distribution . */
void MCMC::Impl::ImplInfecProfile::build(const dict_str_str &deflist)
{
	if (!deflist.empty())
	{
		// check entries
		tools::check_for_unused_entries(
			tools::headers_from_deflist(deflist), headers_);

		// check if definition is generic or augmented
		str_vec types = tools::parse_list_args(deflist.at("arg_at"));
		types += {deflist.at("max_at")};

		bool is_augm = false;
		for (const std::string &at : types)
			if (at == TYPE_AUGM)
				is_augm = true;

		if (is_augm)
			build_augmented_profile_(deflist);
		else
			build_generic_profile_(deflist);
	}
}

/* Build for generic infectivity profile. */
void MCMC::Impl::ImplInfecProfile::build_generic_profile_(const dict_str_str &deflist)
{
	augm_infec_profile_deflist_.clear();
	inst_ = &ImplInfecProfile::do_inst_default_;
	cumul_ = &ImplInfecProfile::do_cumul_default_;

	//----- Name -----//
	if ((deflist.count("name")) && (!tools::is_na(deflist.at("name"))))
	{
		profile_ = compose_distr_registries(
			deflist.at("name"),
			"infectivity profile",
			{&InfecProfileRegistry::inst(),
			 &ContinuosDistrRegistry::inst(),
			 &DiscreteDistrRegistry::inst()});
	}
	else
		throw std::invalid_argument("Missing name.");

	//----- Discrete -----//
	bool is_discrete = true;
	if (!profile_->is_discrete()) // not defined discrete
		if (deflist.count("discr"))
			is_discrete = tools::stobool(deflist.at("discr"));

	//----- Arguments -----//
	if (deflist.count("arg_key") && deflist.count("arg_val") && deflist.count("arg_at"))
		assign_arguments_(tools::parse_list_args(deflist.at("arg_key")),
						  tools::parse_list_args(deflist.at("arg_val")),
						  tools::parse_list_args(deflist.at("arg_at")));

	//----- Range minimum -----//
	assign_single_attribute_("min",
							 (is_discrete) ? "1" : "0", TYPE_SCALAR,
							 HSP_ATTR_RANGE);

	//----- Range maximum -----//
	if (deflist.count("max") && deflist.count("max_at"))
		assign_single_attribute_("max",
								 deflist.at("max"), deflist.at("max_at"),
								 HSP_ATTR_RANGE);

	//----- Reshaping -----//
	profile_->init_reshaping(is_discrete);

	// if discrete and generic, initializes lookup table
	init_lookup_();
}

/* Build for individual augmented function infectivity profile. */
void MCMC::Impl::ImplInfecProfile::build_augmented_profile_(const dict_str_str &deflist)
{
	profile_ = nullptr;
	inst_ = &ImplInfecProfile::do_inst_augm_;
	cumul_ = &ImplInfecProfile::do_cumul_augm_;

	std::string func_name;
	if ((deflist.count("name")) && (!tools::is_na(deflist.at("name"))))
		func_name = deflist.at("name");
	else
		throw std::invalid_argument("Missing name.");

	// auxiliary vector of string with single NA option for generation time
	const str_vec aux_str({HSP_NA_STRING});

	// adjust arguments to connect to generation time
	std::string arg_key = deflist.at("arg_key");
	if (!tools::is_na(arg_key))
		arg_key += C_LIST_DELIM;

	std::string arg_val = deflist.at("arg_val");
	if (!tools::is_na(arg_val))
		arg_val += C_LIST_DELIM;

	std::string arg_at = deflist.at("arg_at");
	if (!tools::is_na(arg_at))
		arg_at += C_LIST_DELIM;

	// range minimum
	if (deflist.count("min") || deflist.count("min_at"))
		throw std::invalid_argument(
			"Impossible to provide custom minimum for infectivity profile.");

	// initialize for instantaneous and cumulative profiles
	bool is_discrete = tools::stobool(deflist.at("discr"));

	for (auto &entry : pImpl_.augm.file_headers)
		augm_infec_profile_deflist_.emplace(entry, str_vec(2, HSP_NA_STRING));

	augm_infec_profile_deflist_.at("const") = str_vec(2, tools::truth(false));
	augm_infec_profile_deflist_.at("discr") = str_vec(2, tools::truth(is_discrete));
	augm_infec_profile_deflist_.at("init_rand") = str_vec(2, tools::truth(false));

	// match infectivity profile and augmented variable definition
	const dict_str_vec matcher = {
		{"name", str_vec({HSP_PROFILE_INST, HSP_PROFILE_CUMUL})},
		{"func", str_vec({"d" + func_name, "p" + func_name})},
		{"farg_key", str_vec(2, arg_key)},
		{"farg_val", str_vec(2, arg_val)},
		{"farg_at", str_vec(2, arg_at)},
		{"min", str_vec(2, (is_discrete) ? "1" : "0")},
		{"min_at", str_vec(2, TYPE_SCALAR)},
		{"max", str_vec(2, deflist.at("max"))},
		{"max_at", str_vec(2, deflist.at("max_at"))}};

	for (auto &pair : matcher)
		augm_infec_profile_deflist_.at(pair.first) = pair.second;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  MCMC : Impl : InfecProfile : getter methods  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Default method for instantaneous infectivity profile. */
double MCMC::Impl::ImplInfecProfile::do_inst_default_(const Individual *, double t)
{
	if (profile_->valid(t))
		return profile_->pdf(t);
	else
		return 0.;
}
/* Default method for cumulative infectivity profile. */
double MCMC::Impl::ImplInfecProfile::do_cumul_default_(const Individual *, double t)
{
	if (profile_->below(t))
		return 0.;
	else if (profile_->above(t))
		return profile_->cdf(profile_->get_max());
	else
		return profile_->cdf(t);
}

/* Lookup method for discrete instantaneous infectivity profile. */
double MCMC::Impl::ImplInfecProfile::do_inst_lookup_(const Individual *, double t)
{
	if (profile_->valid(t) && (t < size_lut_))
		return lookup_inst_.at(static_cast<size_t>(t));
	else
		return 0;
}
/* Lookup method for discrete cumulative infectivity profile. */
double MCMC::Impl::ImplInfecProfile::do_cumul_lookup_(const Individual *, double t)
{
	if (profile_->below(t))
		return 0.;
	else if ((t >= size_lut_) || (profile_->above(t)))
		return lookup_cumul_.at(size_lut_ - 1);
	else
		return lookup_cumul_.at(static_cast<size_t>(t));
}

/* Value of the instantaneous infectivity profile. */
double MCMC::Impl::ImplInfecProfile::do_inst_augm_(const Individual *member, double t)
{
	return member->get_augm_variable(HSP_PROFILE_INST)(t);
}

/* Value of the cumulative infectivity profile. */
double MCMC::Impl::ImplInfecProfile::do_cumul_augm_(const Individual *member, double t)
{
	return member->get_augm_variable(HSP_PROFILE_CUMUL)(t);
}
