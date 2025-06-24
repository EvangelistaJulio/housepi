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
					   indiv hhindex,
					   const std::string &id,
					   int status,
					   const data_key_map<double> &data)
{
	_index = index;
	_hhindex = hhindex;
	_id = id;
	_status = (status < MISSING_STATUS) ? MISSING_STATUS : status;
	_data_map = data;
}

/* Individual deconstructor. */
Individual::~Individual() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Individual : getter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Value of loglikelihood associated to augmented variable models within individual.
 * Store the loglikelihood of individual augmented variables. 
 */
double Individual::compute_and_set_augm_model_loglike()
{
	double x = 0.0;
	for (const std::function<double()> &f : _augm_model_logpdfs)
		x += f();

	if (tools::is_na(x))
		throw std::runtime_error(
			"Individual \"" + _id + "\" augmented loglikelihood return NaN value.");

	// set computed value and store old in buffer
	_augm_buffer = _augm_loglike;
	_augm_loglike = x;

	return x;
}

/* Individual dataset keys. */
const std::vector<data_key> Individual::keys_dataset() const
{
	std::vector<data_key> keys = {};
	for (const auto &pair : _data_map)
		keys.push_back(pair.first);
	return keys;
}
/* Individual datalist keys. */
const std::vector<data_key> Individual::keys_datalist() const
{
	std::vector<data_key> keys = {};
	for (const auto &pair : *_datalist_ptr)
		keys.push_back(pair.first);
	return keys;
}
/* Individual group rate keys. */
const std::vector<rate_key> Individual::keys_group_rate() const
{
	std::vector<data_key> keys = {};
	for (const auto &pair : _rate_map)
		keys.push_back(pair.first);
	return keys;
}
/* Individual augmented variable keys. */
const std::vector<augm_key> Individual::keys_augm() const
{
	std::vector<augm_key> keys = {};
	for (const auto &pair : _augm_map)
		keys.push_back(pair.first);
	return keys;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Individual : setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Individual::assign_all_group_rates(const rate_key_mapvec<const double *> &ptrs)
{
	for (const auto &pair : ptrs)
	{
		rate_key r = pair.first;
		std::vector<const double *> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate \"" + r + "\" for individual " + _id + ".");
		else
			this->_rate_map[r] = GroupRate(r, vec_par_ptr);
	}
}

/* Assign definitions for the augmented variable model. */
void Individual::assign_augm(const augm_key &a, const std::string &name,
							 double value, bool is_discrete, bool is_constant,
							 const std::string &func_name,
							 const std::string &model_name)
{
	Variable &this_augm = _augm_map.at(a);
	this_augm.init(name, value, is_discrete, is_constant, func_name);
	this_augm.assign_distr(model_name);
}

/* Associate functional methods for augmented model pdfs. */
void Individual::init_augm_model(const std::vector<augm_key> &augm_model_list)
{
	_augm_model_logpdfs.clear();
	for (const augm_key &a : augm_model_list)
	{
		const Variable &xaugm = _augm_map.at(a);
		_augm_model_logpdfs.emplace_back([&xaugm]()
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
								 const std::vector<const double *> &pointers)
{
	if (!tools::is_na(name))
	{
		if (pointers.empty())
			throw std::invalid_argument(
				"No parameters associated to group rate \"" + name + "\".");

		this->group_rate_name = name;
		this->n_rates = pointers.size();
		this->pointers = new const double *[pointers.size()];
		for (ushort c = 0; c < this->n_rates; ++c)
			this->pointers[c] = pointers.at(c);
	}
}
Individual::GroupRate::~GroupRate() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Individual::Group Rate : getter  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Individual::GroupRate::get() const
{
	double output = 1.0;
	for (ushort r = 0; r < n_rates; ++r)
		output *= *pointers[r];

	return output;
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
	_house_id = house_id;
	_index = H;
	_house_multiplicity = 1;
	_multiplicity_class = H;
	_house_size = 0;
}
Household::~Household() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : setter  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign a multiplicity to the household. */
void Household::set_multiplicity(ulong m, size_t class_idx)
{
	_house_multiplicity = m;
	_multiplicity_class = class_idx;
}

/* Link an individual to the household. */
void Household::assign_member(size_t index,
							  const std::string &id,
							  int status,
							  const data_key_map<double> &data)
{
	//----- update status lists within the household -----//
	if (status >= INFECTED_STATUS)
		_infec_list.push_back(_house_size);

	else if (status == HEALTHY_STATUS)
		_healthy_list.push_back(_house_size);

	else
		_missing_list.push_back(_house_size);

	//----- emplace new individual -----//
	_members.emplace_back(index, _house_size, id, status, data);
	++_house_size;
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household::change_status(Individual &member, int new_status)
{
	int old_status = member.status();

	if (new_status == old_status)
		return false;

	else
	{
		indiv i = member.hhindex();

		if (old_status >= INFECTED_STATUS)
			_infec_list.erase(std::find(_infec_list.begin(), _infec_list.end(), i));

		else if (old_status == HEALTHY_STATUS)
			_healthy_list.erase(std::find(_healthy_list.begin(), _healthy_list.end(), i));

		else if (old_status == MISSING_STATUS)
			_missing_list.erase(std::find(_missing_list.begin(), _missing_list.end(), i));

		else
			throw std::runtime_error(
				"Unknown status \"" + std::to_string(old_status) + "\" for individual " + member.id() + ".");

		if (new_status >= INFECTED_STATUS)
			_infec_list.push_back(i);

		else if (new_status == HEALTHY_STATUS)
			_healthy_list.push_back(i);

		else if (new_status == MISSING_STATUS)
			_missing_list.push_back(i);

		else
			throw std::runtime_error(
				"Unknown new status \"" + std::to_string(new_status) + "\" for individual " + member.id() + ".");

		member.set_status(new_status);
		return true;
	}
}

/* Associate a pointer to parameter value for the contact from->to.
Allocate memory if it  */
void Household::assign_contact_rate(indiv from, indiv to, const double *x)
{
	// initialize contact rate matrix if necessary
	if (_contact_rates.empty())
		_contact_rates.resize(size(), std::vector<const double *>(size(), &REFERENCE_VALUE));

	_contact_rates[from][to] = x;
}

/* Store the loglikelihood of household observation. */
void Household::set_obs_loglike(double x)
{
	if (tools::is_na(x))
		throw std::runtime_error(
			"Household \"" + _house_id + "\" observation loglikelihood return NaN value.");

	_obs_buffer = _obs_loglike;
	_obs_loglike = x;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Recover previous values of the loglikelihood and cached depedencies. */
void Household::recover_buffer_values(bool is_cache_depen)
{
	_obs_loglike = _obs_buffer;
	if (is_cache_depen)
	{
		// swap buffer into  cache attributes
		_cache_int_vec.swap(_buffer_cache_int_vec);
		_cache_dbl_vec.swap(_buffer_cache_dbl_vec);
		_cache_dbl_mtx.swap(_buffer_cache_dbl_mtx);
	}
}

/* Deep copies of cache attributes in the buffer */
void Household::buffer_all_cache()
{
	_buffer_cache_int_vec = _cache_int_vec;
	_buffer_cache_dbl_vec = _cache_dbl_vec;
	_buffer_cache_dbl_mtx = _cache_dbl_mtx;
}
