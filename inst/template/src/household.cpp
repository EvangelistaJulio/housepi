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
					   const std::string &id,
					   int status,
					   const dataIDmap<double> &data)
{
	_index = index;
	_id = id;
	_status = (status < MISSING_STATUS) ? MISSING_STATUS : status;
	_data_map = data;
}

/* Individual deconstructor. */
Individual::~Individual() = default;

const std::string &Individual::id() const { return _id; }
size_t Individual::index() const { return _index; }
int Individual::status() const { return _status; }

double Individual::data(const dataID &d) const { return _data_map.at(d); }
double Individual::group_rate(const rateID &r) const { return _rate_map.at(r).get(); }
double Individual::augm(const augmID &a) const { return _augm_map.at(a).get(); }
Variable &Individual::get_augm_variable(const augmID &a) { return _augm_map.at(a); }

double Individual::get_augm_model_loglike() const
{
	double indiv_loglike = 0.0;
	for (const std::function<double()> &f : _augm_model_logpdfs)
		indiv_loglike += f();
	return indiv_loglike;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Individual : setter  //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Individual::set_status(int s) { _status = s; }

void Individual::assign_all_group_rates(const rateIDmapvec<const double *> &ptrs)
{
	for (const auto &pair : ptrs)
	{
		rateID r = pair.first;
		std::vector<const double *> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate \"" + r + "\" for individual " + _id + ".");
		else
			this->_rate_map[r] = GroupRate(r, vec_par_ptr);
	}
}

void Individual::clear_augm(const augmID &a) { _augm_map[a] = Variable(); }

/* Assign definitions for the augmented variable model. */
void Individual::assign_augm(const augmID &a, const std::string &name,
							 double value, bool is_discrete, bool is_constant,
							 const std::string &func_name,
							 const std::string &model_name)
{
	Variable &this_augm = _augm_map.at(a);
	this_augm.init(name, value, is_discrete, is_constant, func_name);
	this_augm.assign_distr(model_name);
}

/* Associate functional methods for augmented model pdfs. */
void Individual::init_augm_model(const std::vector<augmID> &augm_model_list)
{
	_augm_model_logpdfs.clear();
	for (const augmID &a : augm_model_list)
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
							  const dataIDmap<double> &data)
{
	//----- update status lists within the household -----//
	if (status >= INFECTED_STATUS)
	{
		++_n_infec;
		_infec_list.push_back(_house_size);
	}
	else if (status == HEALTHY_STATUS)
	{
		++_n_healthy;
		_healthy_list.push_back(_house_size);
	}
	else
	{
		++_n_missing;
		_missing_list.push_back(_house_size);
	}

	//----- emplace new individual -----//
	_members.emplace_back(index, id, status, data);
	++_house_size;
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household::change_status(indiv i, int new_status)
{
	Individual &this_indiv = _members[i];
	int old_status = this_indiv.status();

	if (new_status == old_status)
		return false;

	else
	{
		if (old_status >= INFECTED_STATUS)
		{
			--_n_infec;
			_infec_list.erase(std::find(_infec_list.begin(), _infec_list.end(), i));
		}
		else if (old_status == HEALTHY_STATUS)
		{
			--_n_healthy;
			_healthy_list.erase(std::find(_healthy_list.begin(), _healthy_list.end(), i));
		}
		else if (old_status == MISSING_STATUS)
		{
			--_n_missing;
			_missing_list.erase(std::find(_missing_list.begin(), _missing_list.end(), i));
		}
		else
			throw std::runtime_error(
				"Unknown status \"" + std::to_string(old_status) + "\" for individual " + this_indiv.id() + ".");

		if (new_status >= INFECTED_STATUS)
		{
			++_n_infec;
			_infec_list.push_back(i);
		}
		else if (new_status == HEALTHY_STATUS)
		{
			++_n_healthy;
			_healthy_list.push_back(i);
		}
		else if (new_status == MISSING_STATUS)
		{
			++_n_missing;
			_missing_list.push_back(i);
		}
		else
			throw std::runtime_error(
				"Unknown new status \"" + std::to_string(new_status) + "\" for individual " + this_indiv.id() + ".");

		this_indiv.set_status(new_status);
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
	if ((tools::is_na(_obs_loglike)) && (std::isinf(_obs_loglike)))
		throw std::runtime_error(
			"NaN value found in household observation loglikelihood.");
	if (tools::is_na(x))
		throw std::runtime_error(
			"NaN value set as household observation loglikelihood.");

	_obs_buffer = _obs_loglike;
	_obs_loglike = x;
}

/* Store the loglikelihood of household augmented variables. */
void Household::set_augm_loglike(double x)
{
	if ((tools::is_na(_augm_loglike)) && (std::isinf(_augm_loglike)))
		throw std::runtime_error(
			"NaN value found in household augmented variables loglikelihood.");
	if (tools::is_na(x))
		throw std::runtime_error(
			"NaN value set as household augmented variables loglikelihood.");

	_augm_buffer = _augm_loglike;
	_augm_loglike = x;
}

/* Restore previous values of the loglikelihood. */
void Household::recover_loglike_value()
{
	_obs_loglike = _obs_buffer;
	_augm_loglike = _augm_buffer;
}

//===== GETTER =====//
const std::string &Household::house_id() const { return _house_id; }
size_t Household::house_index() const { return _index; }
ushort Household::size() const { return _house_size; }

const std::string &Household::id(indiv i) const { return _members[i].id(); }
size_t Household::index(indiv i) const { return _members[i].index(); }

ulong Household::multiplicity() const { return _house_multiplicity; }
size_t Household::multiplicity_class_index() const { return _multiplicity_class; }
bool Household::is_multiplicity_head_class() const { return (_index == _multiplicity_class); }

int Household::status(indiv i) const { return _members[i].status(); }
double Household::contact_rate(indiv from, indiv to) const { return *_contact_rates[from][to]; }
Individual &Household::get_member(indiv i) { return _members[i]; }

double Household::data(indiv i, const dataID &d) const { return _members[i].data(d); }
double Household::group_rate(indiv i, const rateID &r) const { return _members[i].group_rate(r); }
double Household::augm(indiv i, const augmID &a) const { return _members[i].augm(a); }

ushort Household::n_infec_status() const { return _n_infec; }
ushort Household::n_healthy_status() const { return _n_healthy; }
ushort Household::n_missing_status() const { return _n_missing; }

double Household::get_augm_loglike() const { return _augm_loglike; }
double Household::get_obs_loglike() const { return _obs_loglike; }

const std::vector<indiv> &Household::infected_status_list() const { return _infec_list; }
const std::vector<indiv> &Household::healthy_status_list() const { return _healthy_list; }
const std::vector<indiv> &Household::missing_status_list() const { return _missing_list; }

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : setter  ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Household::assign_all_group_rates(indiv i, const rateIDmapvec<const double *> &x)
{
	_members[i].assign_all_group_rates(x);
}

/* Set a value for the augmented variable. */
void Household::set_augm_value(indiv i, const augmID &a, double x)
{
	_members[i].get_augm_variable(a).set_value(x);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Household : augmented moves  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Update the values for posterior sampling. */
void Household::init_diagnostics(ulong sample_size)
{
	_smpl_loglike.clear();
	_smpl_loglike = dbl_vec(sample_size);
}

/* Update the values for posterior sampling. */
void Household::update_diagnostics(ulong c)
{
	_smpl_loglike[c] = _obs_loglike / _house_multiplicity + _augm_loglike;
}

/* Return the loglikelihood sampled vector. */
const dbl_vec &Household::get_smpl_loglike() const { return _smpl_loglike; }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : cache attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if any cache quantity is associated to the household. */
bool Household::is_any_cache()
{
	return !(_cache_int_vec.empty() && _cache_dbl_vec.empty() && _cache_dbl_mtx.empty());
}

/* Deep copies of cache attributes in the buffer */
void Household::buffer_all_cache()
{
	_buffer_cache_int_vec = _cache_int_vec;
	_buffer_cache_dbl_vec = _cache_dbl_vec;
	_buffer_cache_dbl_mtx = _cache_dbl_mtx;
}

/* Swap buffer into  cache attributes */
void Household::recover_all_cache()
{
	_cache_int_vec.swap(_buffer_cache_int_vec);
	_cache_dbl_vec.swap(_buffer_cache_dbl_vec);
	_cache_dbl_mtx.swap(_buffer_cache_dbl_mtx);
}

//----- Extra integer vector -----//
void Household::allocate_cache_integer_vector(const cacheID &e) { _cache_int_vec[e] = int_vec(size(), 0); }
int_vec &Household::cache_integer_vector(const cacheID &e) { return _cache_int_vec.at(e); }

//----- Extra double vector -----//
void Household::allocate_cache_double_vector(const cacheID &e) { _cache_dbl_vec[e] = dbl_vec(size(), 0.); }
dbl_vec &Household::cache_double_vector(const cacheID &e) { return _cache_dbl_vec.at(e); }

//----- Extra double matrix -----//
void Household::allocate_cache_double_matrix(const cacheID &e) { _cache_dbl_mtx[e] = matrix<double>(size(), dbl_vec(size(), 0.)); }
matrix<double> &Household::cache_double_matrix(const cacheID &e) { return _cache_dbl_mtx.at(e); }
