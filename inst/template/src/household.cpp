/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * household.cpp
 *
 */

#include <algorithm>

#include "variable.h"
#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Individual : constructors  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Individual constructor. */
Individual::Individual(size_t index,
					   size_t hhindex,
					   const std::string &id,
					   int status,
					   const data_key_map<double> &data)
{
	index_ = index;
	hhindex_ = hhindex;
	id_ = id;
	status_ = (status < MISSING_STATUS) ? MISSING_STATUS : status;
	data_map_ = data;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Individual : getter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Value of loglikelihood associated to augmented variable models within individual.
 * Store the loglikelihood of individual augmented variables.
 */
double Individual::compute_and_set_augm_model_loglike()
{
	double x = 0.0;
	for (const std::function<double()> &f : augm_model_logpdfs_)
		x += f();

	if (tools::is_na(x))
		throw std::runtime_error(
			"Individual \"" + id_ + "\" augmented loglikelihood returned NaN.");

	// set computed value and store old in buffer
	augm_buffer_ = augm_loglike_;
	augm_loglike_ = x;

	return x;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Individual : setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Individual::init_all_group_rates(const rate_key_mapvec<const double *> &ptrs)
{
	for (const auto &pair : ptrs)
	{
		rate_key r = pair.first;
		std::vector<const double *> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate \"" + r + "\" for individual " + id_ + ".");
		else
			this->rate_map_[r] = GroupRate(r, vec_par_ptr);
	}
}

/* Assign definitions for the augmented variable model. */
void Individual::init_augm(const augm_key &a, const std::string &name,
							 double value, bool is_discrete, bool is_constant,
							 const std::string &func_name,
							 const std::string &model_name)
{
	Variable &this_augm = safe_at(augm_map_, a);
	this_augm.init(name, value, is_discrete, is_constant, func_name);
	this_augm.init_distr(model_name);
}

/* Associate functional methods for augmented model pdfs. */
void Individual::init_augm_model(const std::vector<augm_key> &augm_model_list)
{
	augm_model_logpdfs_.clear();
	for (const augm_key &a : augm_model_list)
	{
		const Variable &xaugm = safe_at(augm_map_, a);
		augm_model_logpdfs_.emplace_back([&xaugm]()
										 { return xaugm.distr_logpdf(); });
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Individual::Group Rate : constructor  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

Individual::GroupRate::GroupRate(const std::string &name,
								 const std::vector<const double *> &ptrs)
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

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Household : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

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

/* Assign a multiplicity to the household. */
void Household::set_multiplicity(ulong m, size_t class_idx)
{
	house_multiplicity_ = m;
	multiplicity_class_ = class_idx;
}

/* Link an individual to the household. */
void Household::init_member(size_t index,
							  const std::string &id,
							  int status,
							  const data_key_map<double> &data)
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

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
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

/* Associate a pointer to parameter value for the contact from->to.
Allocate memory if it  */
void Household::init_contact_rates(size_t from, size_t to, const double *x)
{
	// initialize contact rate matrix if necessary
	if (contact_rates_.empty())
		contact_rates_.resize(size(), std::vector<const double *>(size(), &REFERENCE_VALUE));

	contact_rates_[from][to] = x;
}

/* Store the loglikelihood of household observation. */
void Household::set_obs_loglike(double x)
{
	if (tools::is_na(x))
		throw std::runtime_error(
			"Household \"" + house_id_ + "\" observation loglikelihood returned NaN.");

	obs_buffer_ = obs_loglike_;
	obs_loglike_ = x;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Recover previous values of the loglikelihood and cached depedencies. */
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

/* Deep copies of cache attributes in the buffer */
void Household::buffer_all_cache()
{
	buffer_cache_int_vec_ = cache_int_vec_;
	buffer_cache_size_t_vec_ = cache_size_t_vec_;
	buffer_cache_dbl_vec_ = cache_dbl_vec_;
	buffer_cache_dbl_mtx_ = cache_dbl_mtx_;
}
