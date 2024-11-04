/*
 * household.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

//===== CONSTRUCTOR =====//
Household::Household() {};
Household::Household(const std::string &ID, int H)
{
	this->ID = ID;
	this->index = H;
}
Household::~Household() = default;

/* Link an individual to the household */
void Household::assign_individual(const std::string &ID, int status, const dictargs &data)
{
	//----- update status lists within the household -----//
	if (status >= INFECTIOUS)
	{
		++num_infectious;
		infectious_set.insert(house_size);

		_is_infectious_.push_back(true);
		_is_susceptible_.push_back(false);
		_is_missing_status_.push_back(false);
	}
	else if (status == SUSCEPTIBLE)
	{
		++num_susceptible;
		susceptible_set.insert(house_size);

		_is_infectious_.push_back(false);
		_is_susceptible_.push_back(true);
		_is_missing_status_.push_back(false);
	}
	else if (status == MISSING_STATUS)
	{
		++num_missing_status;
		missing_status_set.insert(house_size);

		_is_infectious_.push_back(false);
		_is_susceptible_.push_back(false);
		_is_missing_status_.push_back(true);
	}
	else
	{
		throw std::invalid_argument(
			"Unknown status '" + std::to_string(status) + "' for individual " + ID + ".");
	}

	//----- update household size -----//
	++house_size;

	//----- emplace new individual -----//
	indiv.emplace_back(ID, status, data);
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household::change_status(int i, int new_status)
{
	Individual &this_indiv = indiv[i];
	int old_status = this_indiv.get_status();

	if (new_status == old_status)
	{
		return false;
	}
	else
	{
		if (old_status >= INFECTIOUS)
		{
			--num_infectious;
			infectious_set.erase(i);
			_is_infectious_[i] = false;
		}
		else if (old_status == SUSCEPTIBLE)
		{
			--num_susceptible;
			susceptible_set.erase(i);
			_is_susceptible_[i] = false;
		}
		else if (old_status == MISSING_STATUS)
		{
			--num_missing_status;
			missing_status_set.erase(i);
			_is_missing_status_[i] = false;
		}
		else
		{
			std::ostringstream oss;
			oss << "While changing status: unknown old status '"
				<< old_status << "' of individual " << this_indiv.get_ID() << ".";
			throw std::runtime_error(oss.str());
		}

		if (new_status >= INFECTIOUS)
		{
			++num_infectious;
			infectious_set.insert(i);
			_is_infectious_[i] = true;
		}
		else if (new_status == SUSCEPTIBLE)
		{
			++num_susceptible;
			susceptible_set.insert(i);
			_is_susceptible_[i] = true;
		}
		else if (new_status == MISSING_STATUS)
		{
			++num_missing_status;
			missing_status_set.insert(i);
			_is_missing_status_[i] = true;
		}
		else
		{
			std::ostringstream oss;
			oss << "While changing status: unknown new status '"
				<< new_status << "' for individual " << this_indiv.get_ID() << ".";
			throw std::runtime_error(oss.str());
		}

		this_indiv.set_status(new_status);
		return true;
	}
}

/* Initializes matrix of contact rates. */
void Household::init_contact_matrix()
{
	this->_which_rate_.resize(size(), std::vector<std::string>(size(), UNDEFINED));
	this->contact_rates.resize(size(), std::vector<const double *>(size(), &REFERENCE_VALUE));
}

void Household::assign_contact_rate(int from, int to, const std::string &parID) { this->_which_rate_[from][to] = parID; }
void Household::set_contact_rate(int from, int to, const double *x) { this->contact_rates[from][to] = x; }
void Household::set_loglike(double x)
{
	this->buffer = this->house_loglike;
	this->house_loglike = x;
}
void Household::restore_loglike() { this->house_loglike = this->buffer; }

//===== GETTER =====//
const std::string Household::get_ID() const { return this->ID; }
const std::string Household::get_ID(int i) const { return this->indiv[i].get_ID(); }
int Household::get_index() const { return this->index; }
int Household::size() const { return this->house_size; }

std::unordered_map<std::string, std::string> Household::get_group_rate_map(int i) const { return indiv[i].get_group_rate_map(); }
double Household::get_group_rate(int i, const std::string &rateID) const { return this->indiv[i].get_group_rate(rateID); }

const std::string Household::get_contact_rate_map(int from, int to) const { return this->_which_rate_[from][to]; }
double Household::get_contact_rate(int from, int to) const { return *this->contact_rates[from][to]; }

int Household::get_status(int i) const { return this->indiv[i].get_status(); }
double Household::get_augm(int i, const std::string &augmID) const { return this->indiv[i].get_augm(augmID); }
double Household::get_data(int i, const std::string &dataID) const { return this->indiv[i].get_data(dataID); }
const double *Household::get_data_pointer(int i, const std::string &dataID) const { return this->indiv[i].get_data_pointer(dataID); }

int Household::get_num_infectious() const { return this->num_infectious; }
int Household::get_num_susceptible() const { return this->num_susceptible; }
int Household::get_num_missing_status() const { return this->num_missing_status; }

const std::set<int> &Household::get_infectious_set() const { return this->infectious_set; }
const std::set<int> &Household::get_susceptible_set() const { return this->susceptible_set; }
const std::set<int> &Household::get_missing_status_set() const { return this->missing_status_set; }

double Household::get_loglike() const { return this->house_loglike; }
double Household::get_delta_loglike() const { return this->house_loglike - this->buffer; }

double Household::augm_model_logpdf(int i, const std::string &augmID) const { return this->indiv[i].augm_model_logpdf(augmID); }
double Household::augm_model_logratio(int i, const std::string &augmID) const { return this->indiv[i].augm_model_logratio(augmID); }
double Household::augm_propose_logratio(int i, const std::string &augmID) const { return this->indiv[i].augm_propose_logratio(augmID); }
double Household::augm_move_logratio(int i, const std::string &augmID) const { return this->indiv[i].augm_move_logratio(augmID); }

//===== SETTER =====//
void Household::disable_assert()
{
	for (int i = 0; i < size(); ++i)
	{
		indiv[i].disable_assert();
	}
}

void Household::assign_augm(int i, const std::string &augmID, double value,
							const std::string &model_name, dictargs_ptr model_args_ptr, dictargs range,
							const std::string &rule_name, dictargs rule_args, double rule_from)
{
	this->indiv[i].assign_augm(augmID, value,
							   model_name, model_args_ptr, range,
							   rule_name, rule_args, rule_from);
}
void Household::set_augm_value(int i, const std::string &augmID, double x) { this->indiv[i].set_augm_value(augmID, x); }
void Household::augm_model_rand(int i, const std::string &augmID, sys_rnd_gen &prng) { this->indiv[i].augm_model_rand(augmID, prng); }
bool Household::propose_valid_augm_move(int i, const std::string &augmID, sys_rnd_gen &prng) { return this->indiv[i].propose_valid_augm_move(augmID, prng); }
void Household::augm_restore_value(int i, const std::string &augmID) { this->indiv[i].augm_restore_value(augmID); }

void Household::assign_group_rate(int i, const std::string &rateID, const std::string &parID) { this->indiv[i].assign_group_rate(rateID, parID); }
void Household::set_rate(int i, const std::string &rateID, const double *ptr) { this->indiv[i].set_rate(rateID, ptr); }
void Household::clear_augm(int i, const std::string &augmID) { this->indiv[i].clear_augm(augmID); }

//===== DIAGNOSTICS =====//

/* update the values for the Bayesian diagnostics */
void Household::init_diagnostics(ulong sample_size)
{
	dgn_count = 0;
	dgn_loglike.clear();
	dgn_loglike = std::vector<double>(sample_size);
}

/* update the values for the Bayesian diagnostics */
void Household::update_diagnostics()
{
	dgn_loglike[dgn_count++] = house_loglike;
}