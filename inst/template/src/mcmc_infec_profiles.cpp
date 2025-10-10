/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_infec_profiles.cpp
 *
 */

#include "tools/params.h"
#include "tools/read.h"
#include "tools/string.h"
#include "mcmc_impl.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : Impl : InfecProfile : assign template  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Set model template while constructiing model.
void MCMC::Impl::InfecProfile::init_template(u_dict<std::string> deflist)
{
	// if not previous build, proceed with template
	// if no template, assume no need for an infectivity profile
	if (hsp::na::is(profile_->name()) && augm_infec_profile_deflist_.empty())
		if (!deflist.empty())
			try
			{
				if (!(deflist % "func"))
					throw std::invalid_argument(
						"Missing function name for infectivity profile template.");

				const str_vec file_headers = {"func", "discr", hsp::params::MAX};

				for (const std::string &key : file_headers)
					if (!(deflist % key))
						deflist.emplace(key, hsp::na::STRING);

				std::string str_is_discr = "false";
				if (deflist % "discr")
					if (!hsp::na::is(deflist.at("discr")))
						str_is_discr = deflist.at("discr");
				deflist["discr"] = str_is_discr;

				build(deflist);

				// forbid augmented infectivity profile for template
				if (is_augmented())
					throw std::invalid_argument(
						"Infectivity profile templates depending on augmented variables are not supported.");

				// UPDATE: maybe also not depending on parameters?
			}
			catch (const std::exception &e)
			{
				throw std::runtime_error(
					"While building infectivity profile template:\n" + std::string(e.what()));
			}
}

/// Append infectivity profile to augmented variable definitions.
/// UPDATE: problematic in case of redefinitions of infectivity profile
void MCMC::Impl::InfecProfile::append_augm_definition()
{
	if (is_augmented())
	{
		// check for reserved augmented names
		for (const std::string &a : pImpl_.augm.deftable.at("func"))
			for (const std::string &forbid : {hsp::infprof::INSTANTANEOUS, hsp::infprof::CUMULATIVE})
				if (a == forbid)
					throw std::invalid_argument(
						"Augmented variable \"" + forbid + "\" is reserved. Please, pick another name.");

		// append new definitions
		for (const auto &pair : augm_infec_profile_deflist_)
			pImpl_.augm.deftable.at(pair.first) += pair.second;
	}
}

/// Store values of the infectivity profile.
void MCMC::Impl::InfecProfile::init_lookup_()
{
	inst_ = &InfecProfile::do_inst_lookup_;
	cumul_ = &InfecProfile::do_cumul_lookup_;

	if (hsp::na::is(profile_->get_max()))
		size_lut_ = hsp::infprof::LOOKUP_SIZE;
	else
		size_lut_ = profile_->get_max() + 1;

	lookup_inst_.resize(size_lut_);
	lookup_cumul_.resize(size_lut_);
	lookup_buffer_inst_.resize(size_lut_);
	lookup_buffer_cumul_.resize(size_lut_);

	update(true);
}

/// Clear infectivity profile lookup table.
void MCMC::Impl::InfecProfile::clear_lookup_()
{
	size_lut_ = 0;

	lookup_inst_.clear();
	lookup_inst_.shrink_to_fit();
	lookup_cumul_.clear();
	lookup_cumul_.shrink_to_fit();
	lookup_buffer_inst_.clear();
	lookup_buffer_inst_.shrink_to_fit();
	lookup_buffer_cumul_.clear();
	lookup_buffer_cumul_.shrink_to_fit();
}

/// If infectivity profile is a lookup table, store its values.
void MCMC::Impl::InfecProfile::update(bool is_cache_depen)
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

/// If infectivity profile is a lookup table, recover its buffered values.
void MCMC::Impl::InfecProfile::recover(bool is_cache_depen)
{
	// UPDATE: this should not be recoverd if constant
	if (is_cache_depen)
		if (is_lookup())
		{
			lookup_inst_.swap(lookup_buffer_inst_);
			lookup_cumul_.swap(lookup_buffer_cumul_);
		}
}

/// Build infectivity profile distribution .
void MCMC::Impl::InfecProfile::build(const u_dict<std::string> &deflist)
{
	if (!deflist.empty())
	{
		const str_vec file_headers = {"func", "discr", hsp::params::MAX};

		// check entries
		hsp::str::assert_no_unused_entries(hsp::read::list_headers(deflist), file_headers);

		// check definition list has a function
		if (!(deflist % "func"))
			throw std::invalid_argument(
				"Missing function definition for infectivity profile.");
		else if (DefinitionParser(deflist.at("func")).is_na())
			throw std::invalid_argument(
				"Empty function definition for infectivity profile.");

		// check if definition is generic or augmented
		bool is_augm = DefinitionParser(deflist.at("func")).is_augmented_depen() ||
					   DefinitionParser(deflist.at(hsp::params::MAX)).is_augmented_depen();

		if (is_augm)
			build_augmented_profile_(deflist);
		else
			build_generic_profile_(deflist);
	}
}

/// Build for generic infectivity profile.
void MCMC::Impl::InfecProfile::build_generic_profile_(const u_dict<std::string> &deflist)
{
	augm_infec_profile_deflist_.clear();
	inst_ = &InfecProfile::do_inst_default_;
	cumul_ = &InfecProfile::do_cumul_default_;

	DefinitionParser parser(DefinitionParser(deflist.at("func")));
	profile_ = compose_distr_registries(parser.root.func,
										"infectivity profile",
										{&InfecProfileRegistry::inst(),
										 &ContinuosDistrRegistry::inst(),
										 &DiscreteDistrRegistry::inst()});
	for (const auto &node : parser.root.children)
	{
		profile_->init_attr(node.key, pImpl_.read_attribute_(node));
		pImpl_.make_param_cache_depen_(node);
	}

	// discretization
	if (deflist % "discr")
		if (hsp::str::tobool(deflist.at("discr")))
			profile_->init_discretization();

	// range minimum
	if (profile_->is_discrete())
		profile_->init_attr(hsp::params::MIN, pImpl_.read_attribute_("1|N"));
	else
		profile_->init_attr(hsp::params::MIN, pImpl_.read_attribute_("0|N"));

	// range maximum
	if (deflist % hsp::params::MAX)
		profile_->init_attr(hsp::params::MAX, pImpl_.read_attribute_(deflist.at(hsp::params::MAX)));
	// range truncation
	profile_->init_truncation();

	// if discrete and generic, initializes lookup table
	if (profile_->is_discrete())
		init_lookup_();
	else
		clear_lookup_();
}

/// Build for individual augmented function infectivity profile.
void MCMC::Impl::InfecProfile::build_augmented_profile_(const u_dict<std::string> &deflist)
{
	// use generic build to check profile name at this stage
	build_generic_profile_({{"func", DefinitionParser(deflist.at("func")).root.func}});

	// drop generic definition
	clear_lookup_();
	profile_ = std::make_unique<Distribution>();
	inst_ = &InfecProfile::do_inst_augm_;
	cumul_ = &InfecProfile::do_cumul_augm_;

	// initialize for instantaneous and cumulative profiles
	bool is_discrete = hsp::str::tobool(deflist.at("discr"));

	// match infectivity profile and augmented variable definition
	augm_infec_profile_deflist_ = {
		{"name", str_vec({hsp::infprof::INSTANTANEOUS, hsp::infprof::CUMULATIVE})},
		{"status", str_vec(2, hsp::na::STRING)}, // UPDATE
		{"func", str_vec(2, deflist.at("func"))},
		{"fcumul", {"false", "true"}},
		{"value", str_vec(2, hsp::na::STRING)},
		{"const", str_vec(2, "false")},
		{"discr", str_vec(2, hsp::str::format(is_discrete))},
		{"model", str_vec(2, hsp::na::STRING)},
		{hsp::params::MIN, str_vec(2, (is_discrete) ? "1|N" : "0|N")},
		{hsp::params::MAX, str_vec(2, deflist.at(hsp::params::MAX))},
		{"rule", str_vec(2, hsp::na::STRING)},
		{"scale", str_vec(2, hsp::na::STRING)}};
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  MCMC : Impl : InfecProfile : getter methods  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Lookup method for discrete instantaneous infectivity profile.
double MCMC::Impl::InfecProfile::do_inst_lookup_(const Individual *, double t)
{
	if (profile_->is_valid(t) && (t < size_lut_))
		return lookup_inst_.at(static_cast<size_t>(t));
	else
		return 0;
}
/// Lookup method for discrete cumulative infectivity profile.
double MCMC::Impl::InfecProfile::do_cumul_lookup_(const Individual *, double t)
{
	if (profile_->is_below(t))
		return 0.;
	else if ((t >= size_lut_) || (profile_->is_above(t)))
		return lookup_cumul_.at(size_lut_ - 1);
	else
		return lookup_cumul_.at(static_cast<size_t>(t));
}
