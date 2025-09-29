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
 * household.cpp
 *
 */

#include <algorithm>

#include "variable.h"
#include "household.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  AugmentedMap : setter  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign definitions for the augmented variable model.
void AugmentedMap::init_augm(const std::string &a,
							 double value,
							 bool is_const,
							 bool is_discr,
							 const DefinitionParser &func_name,
							 bool is_prim,
							 const DefinitionParser &model_name)
{
	Variable &this_augm = safe_at(augm_map_, a, "augmented variable");
	this_augm.init(a, value, is_const, is_discr, func_name, is_prim);
	this_augm.init_distr(model_name);
}

/// Associate functional methods for augmented model pdfs.
void AugmentedMap::init_augm_model_logratios(const u_dict<str_vec> &models_ordered)
{
	// >> clean pre-existing
	for (auto &pair : model_logratios_)
	{
		pair.second.clear();
		pair.second.shrink_to_fit();
	}
	model_logratios_.clear();
	model_logratios_.rehash(0);
	// <<

	// associate models if applicable
	for (const auto &pair : models_ordered)
	{
		const std::string &a = pair.first;
		model_logratios_[a] = {};

		// repeat for all dependent models
		for (const std::string &b : pair.second)
		{
			const Variable &xaugm = augm_map_[b];
			model_logratios_[a].emplace_back([&xaugm]()
											 { return xaugm.get_distr_logratio(); });
		}
	}
}

/// Value of loglikelihood associated to augmented variable models.
/// Store the loglikelihood of individual augmented variables.
double AugmentedMap::compute_and_set_augm_model_loglike()
{
	buffer_ = loglike_; // store old value in buffer
	loglike_ = 0.0;

	for (const auto &pair : model_logratios_)
		loglike_ += augm_map_[pair.first].get_distr_logpdf();

	if (tools::is_na(loglike_))
		throw std::runtime_error(
			"Augmented loglikelihood returned NaN for " + err_name() + ".");

	return loglike_;
}

/// Change in loglikelihood associated to given augmented.
double AugmentedMap::compute_and_set_augm_model_logratio(const std::string &a)
{
	buffer_ = loglike_; // store old value in buffer

	double delta_ll = 0.;
	for (const std::function<double()> &f : model_logratios_[a])
		delta_ll += f();
	loglike_ += delta_ll;

	if (tools::is_na(loglike_))
		throw std::runtime_error(
			"Augmented loglikelihood returned NaN for " + err_name() +
			" after changing " + augm_map_[a].err_name() + ".");

	return delta_ll;
}

/// Restore user values if applicable.
void AugmentedMap::user_reset()
{
	for (auto &pair : augm_map_)
		pair.second.user_reset();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Individual : constructor  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Individual constructor.
Individual::Individual(size_t index,
					   size_t hhindex,
					   const std::string &id,
					   int status,
					   const u_dict<double> &data)
{
	index_ = index;
	hhindex_ = hhindex;
	id_ = id;
	status_ = status;
	data_map_ = data;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Individual : setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize all group rates for the individual.
void Individual::init_all_group_rates(const u_dict<std::vector<observ_dbl_ptr>> &ptrs)
{
	for (const auto &pair : ptrs)
	{
		std::string r = pair.first;
		std::vector<observ_dbl_ptr> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate \"" + r + "\" for " + err_name() + ".");
		else
			rate_map_[r] = GroupRate_(r, vec_par_ptr);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Individual::Group Rate : constructor  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Individual::GroupRate_::GroupRate_(const std::string &name,
								   const std::vector<observ_dbl_ptr> &ptrs)
{
	if (!tools::is_na(name))
	{
		if (ptrs.empty())
			throw std::invalid_argument(
				"No parameters associated to group rate \"" + name + "\".");
		name_ = name;
		pointers_ = ptrs;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Household : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Household constructor.
Household::Household(const std::string &house_id, size_t H)
{
	house_id_ = house_id;
	index_ = H;
	house_multiplicity_ = 1;
	multiplicity_class_ = H;
	house_size_ = 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : setter  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Assign a multiplicity to the household.
void Household::set_multiplicity(unsigned long m, size_t class_idx)
{
	house_multiplicity_ = m;
	multiplicity_class_ = class_idx;
}

/// Link an individual to the household.
void Household::init_member(size_t index,
							const std::string &id,
							int status,
							const u_dict<double> &data)
{
	// emplace new individual
	members_.emplace_back(index, house_size_, id, status, data);

	// update status lists within the household
	if (members_.back().is_infected_status())
		infected_list_.push_back(house_size_);
	else if (members_.back().is_healthy_status())
		healthy_list_.push_back(house_size_);
	else
		missing_list_.push_back(house_size_);

	++house_size_;
}

/// Change_status of individual i to new status (without controls).
///	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference.
bool Household::change_status(Individual &member, int new_status)
{
	if (new_status != member.status())
	{
		const size_t i = static_cast<size_t>(member); // get within house index

		if (member.is_infected_status())
			infected_list_.erase(std::find(infected_list_.begin(), infected_list_.end(), i));
		else if (member.is_healthy_status())
			healthy_list_.erase(std::find(healthy_list_.begin(), healthy_list_.end(), i));
		else
			missing_list_.erase(std::find(missing_list_.begin(), missing_list_.end(), i));

		// change member status
		member.set_status(new_status);

		if (member.is_infected_status())
			infected_list_.push_back(i);
		else if (member.is_healthy_status())
			healthy_list_.push_back(i);
		else
			missing_list_.push_back(i);

		return true;
	}
	return false;
}

/// Associate a pointer to parameter value for the contact from->to.
void Household::init_contact_rates(size_t from, size_t to, observ_dbl_ptr x)
{
	// initialize contact rate matrix if necessary
	if (contact_rates_.empty())
		contact_rates_.resize(size(), std::vector<observ_dbl_ptr>(size(), &HSP_REFER_VALUE));

	contact_rates_[from][to] = x;
}

/// Store the loglikelihood of household observation.
void Household::set_obs_loglike(double x)
{
	if (tools::is_na(x))
		throw std::runtime_error(
			"Observation loglikelihood returned NaN " + err_name() + ".");

	obs_buffer_ = obs_loglike_;
	obs_loglike_ = x;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Recover previous values of the loglikelihood and cached depedencies.
void Household::recover_buffer_values(bool is_cache_depen)
{
	obs_loglike_ = obs_buffer_;
	if (is_cache_depen)
	{
		// swap buffer into  cache attributes
		cache_int_vec_.swap(buffer_cache_int_vec_);
		cache_size_t_vec_.swap(buffer_cache_size_t_vec_);
		cache_dbl_vec_.swap(buffer_cache_dbl_vec_);
		cache_dbl_mtx_.swap(buffer_cache_dbl_mtx_);
	}
}

/// Deep copies of cache attributes in the buffer
void Household::buffer_all_cache()
{
	buffer_cache_int_vec_ = cache_int_vec_;
	buffer_cache_size_t_vec_ = cache_size_t_vec_;
	buffer_cache_dbl_vec_ = cache_dbl_vec_;
	buffer_cache_dbl_mtx_ = cache_dbl_mtx_;
}
