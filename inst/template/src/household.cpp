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
Individual::Individual(const std::string &id,
					   int status,
					   const dataIDmap<double> &data)
{
	// associate the new id
	this->id = id;
	// associate the new status
	this->status = (status < MISSING_STATUS) ? MISSING_STATUS : status;
	// associate the new data
	this->data_map = data;
}

/* Individual deconstructor. */
Individual::~Individual() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////  Individual : assign all group rates  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Individual::assign_all_group_rates(const rateIDmapvec<const double *> &ptrs)
{
	for (const auto &pair : ptrs)
	{
		rateID r = pair.first;
		std::vector<const double *> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate \"" + r + "\" for individual " + id + ".");
		else
			this->rate_map[r] = GroupRate(r, vec_par_ptr);
	}
}

double Individual::data(const dataID &d) const { return data_map.at(d); }
double Individual::group_rate(const rateID &r) const { return rate_map.at(r).get(); }
double Individual::augm(const augmID &a) const { return augm_map.at(a).get(); }

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
void Household::assign_individual(const std::string &id,
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
	_members.emplace_back(id, status, data);
	++_house_size;
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household::change_status(indiv i, int new_status)
{
	Individual &this_indiv = _members[i];
	int old_status = this_indiv.status;

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
				"Unknown status \"" + std::to_string(old_status) + "\" for individual " + this_indiv.id + ".");

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
				"Unknown new status \"" + std::to_string(new_status) + "\" for individual " + this_indiv.id + ".");

		this_indiv.status = new_status;
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
	if (std::isnan(x))
		throw std::runtime_error("Set a NaN value as household observation loglikelihood.");

	_obs_buffer = _obs_loglike;
	_obs_loglike = x;
}

/* Store the loglikelihood of household augmented variables. */
void Household::set_augm_loglike(double x)
{
	if (std::isnan(x))
		throw std::runtime_error("Set a NaN value as household augmented variables loglikelihood.");

	_augm_buffer = _augm_loglike;
	_augm_loglike = x;
}

/* Add a local variation to loglikelihood of household augmented variables. */
void Household::increment_augm_loglike(double delta)
{
	set_augm_loglike(_augm_loglike + delta);
}

/* Restore previous values of the loglikelihood. */
void Household::restore_loglike()
{
	_obs_loglike = _obs_buffer;
	_augm_loglike = _augm_buffer;
}

//===== GETTER =====//
const std::string &Household::house_id() const { return _house_id; }
const std::string &Household::id(indiv i) const { return _members[i].id; }
size_t Household::index() const { return _index; }
ushort Household::size() const { return _house_size; }

ulong Household::multiplicity() const { return _house_multiplicity; }
size_t Household::multiplicity_class_index() const { return _multiplicity_class; }
bool Household::is_multiplicity_head_class() const { return (_index == _multiplicity_class); }

int Household::status(indiv i) const { return _members[i].status; }
double Household::contact_rate(indiv from, indiv to) const { return *_contact_rates[from][to]; }

double Household::data(indiv i, const dataID &d) const { return _members[i].data(d); }
double Household::group_rate(indiv i, const rateID &r) const { return _members[i].group_rate(r); }
double Household::augm(indiv i, const augmID &a) const { return _members[i].augm(a); }

double Household::augm_model_logpdf(indiv i, const augmID &a) const { return _members[i].augm_map.at(a).distr_logpdf(); }
double Household::augm_model_logratio(indiv i, const augmID &a) const { return _members[i].augm_map.at(a).distr_logratio(); }
double Household::augm_propose_logratio(indiv i, const augmID &a) const { return _members[i].augm_map.at(a).propose_logratio(); }
double Household::augm_move_logratio(indiv i, const augmID &a) const { return _members[i].augm_map.at(a).get_move_logratio(); }

ushort Household::n_infec_status() const { return _n_infec; }
ushort Household::n_healthy_status() const { return _n_healthy; }
ushort Household::n_missing_status() const { return _n_missing; }

double Household::get_augm_loglike() const { return _augm_loglike; }
double Household::get_obs_loglike() const { return _obs_loglike; }
double Household::get_delta_obs_loglike() const { return _obs_loglike - _obs_buffer; }

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

/* Assign definitions for the augmented variable model. */
void Household::assign_augm(indiv i, const augmID &a, const std::string &name, double value,
							bool is_discrete, const std::string &model_name)
{
	Variable &this_augm = _members[i].augm_map.at(a);
	this_augm.init(name, value, is_discrete);
	this_augm.assign_distribution(model_name);
}

/* Assign scalar argument for augmented variable model. */
void Household::assign_augm_model_scalar_arg(indiv i, const augmID &a, const std::string &arg_name, double x)
{
	_members[i].augm_map.at(a).assign_distr_scalar_arg(arg_name, x);
}
/* Assign pointer to argument for augmented variable model. */
void Household::assign_augm_model_pointer_arg(indiv i, const augmID &a, const std::string &arg_name, const double *y)
{
	_members[i].augm_map.at(a).assign_distr_pointer_arg(arg_name, y);
}

/* Assign scalar range minimum for augmented variable model. */
void Household::assign_augm_model_scalar_range_min(indiv i, const augmID &a, double x)
{
	_members[i].augm_map.at(a).assign_distr_scalar_range_min(x);
}
/* Assign pointer range minimum for augmented variable model. */
void Household::assign_augm_model_pointer_range_min(indiv i, const augmID &a, const double *y)
{
	_members[i].augm_map.at(a).assign_distr_pointer_range_min(y);
}

/* Assign scalar range maximum for augmented variable model. */
void Household::assign_augm_model_scalar_range_max(indiv i, const augmID &a, double x)
{
	_members[i].augm_map.at(a).assign_distr_scalar_range_max(x);
}
/* Assign pointer range maximum for augmented variable model. */
void Household::assign_augm_model_pointer_range_max(indiv i, const augmID &a, const double *y)
{
	_members[i].augm_map.at(a).assign_distr_pointer_range_max(y);
}

void Household::init_augm_model_reshaping(indiv i, const augmID &a)
{
	_members[i].augm_map.at(a).init_distr_reshaping();
}

/* Set a value for the augmented variable. */
void Household::set_augm_value(indiv i, const augmID &a, double x)
{
	_members[i].augm_map.at(a).set_value(x);
}

/* Draw a value from the augmented variable distribution. */
void Household::augm_model_rand(indiv i, const augmID &a, rnd_gen &prng)
{
	_members[i].augm_map.at(a).rand(prng);
}

/* Assign definitions for the augmented variable proposer. */
void Household::assign_augm_proposer(indiv i, const augmID &a, const std::string &rule_name,
									 double rule_scale, double rule_from)
{
	Variable &this_augm = _members[i].augm_map.at(a);
	this_augm.assign_proposer(rule_name, rule_scale, rule_from);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Household : augmented moves  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* True if the propose augm move is valid. */
bool Household::propose_valid_augm_move(indiv i, const augmID &a, rnd_gen &prng)
{
	Variable &this_augm = _members[i].augm_map.at(a);
	this_augm.propose(prng); // propose new value and buffer old

	if (this_augm.valid_move()) // proceed only if valid
		return true;
	else
		this_augm.restore_value(); // reassign buffered value

	return false;
}

/* Restore buffer valies of augmented variable, extra and loglikelihood. */
void Household::restore_augm_value(indiv i, const augmID &a)
{
	_members[i].augm_map.at(a).restore_value();
	restore_all_extra();
	restore_loglike();
}

/* Assign empty augmented variable. */
void Household::clear_augm(indiv i, const augmID &a) { _members[i].augm_map[a] = Variable(); }

/* Update the values for posterior sampling. */
void Household::init_diagnostics(ulong sample_size)
{
	_smpl_loglike.clear();
	_smpl_loglike = dbl_vec(sample_size);
}

/* Update the values for posterior sampling. */
void Household::update_diagnostics(ulong count)
{
	_smpl_loglike[count] = _obs_loglike + _augm_loglike;
}

/* Return the loglikelihood sampled vector. */
const dbl_vec &Household::get_smpl_loglike() const { return _smpl_loglike; }

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Household : extra attributes  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

bool Household::is_any_extra()
{
	return !(_extra_int_vec.empty() && _extra_dbl_vec.empty() && _extra_dbl_mtx.empty());
}

/* Deep copies of extra attributes in the buffer */
void Household::buffer_all_extra()
{
	_buffer_extra_int_vec = _extra_int_vec;
	_buffer_extra_dbl_vec = _extra_dbl_vec;
	_buffer_extra_dbl_mtx = _extra_dbl_mtx;
}

/* Swap buffer into  extra attributes */
void Household::restore_all_extra()
{
	_extra_int_vec.swap(_buffer_extra_int_vec);
	_extra_dbl_vec.swap(_buffer_extra_dbl_vec);
	_extra_dbl_mtx.swap(_buffer_extra_dbl_mtx);
}

//----- Extra integer vector -----//
void Household::allocate_extra_integer_vector(const extraID &e) { _extra_int_vec[e] = int_vec(size(), 0); }
int_vec &Household::extra_integer_vector(const extraID &e) { return _extra_int_vec.at(e); }

//----- Extra double vector -----//
void Household::allocate_extra_double_vector(const extraID &e) { _extra_dbl_vec[e] = dbl_vec(size(), 0.); }
dbl_vec &Household::extra_double_vector(const extraID &e) { return _extra_dbl_vec.at(e); }

//----- Extra double matrix -----//
void Household::allocate_extra_double_matrix(const extraID &e) { _extra_dbl_mtx[e] = matrix<double>(size(), dbl_vec(size(), 0.)); }
matrix<double> &Household::extra_double_matrix(const extraID &e) { return _extra_dbl_mtx.at(e); }
