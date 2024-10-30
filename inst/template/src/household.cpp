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

/* Set all infectious statuses to standard. */
void Household::standardize_infectious_status()
{
	for (int i: infectious_set)
	{
		indiv[i].set_status(INFECTIOUS);
	}
}

/* Initializes matrix of contact rates. */
void Household::init_contact_matrix()
{
	this->_which_rate_.resize(size(), std::vector<std::string>(size(), UNDEFINED));
	this->contact_rates.resize(size(), std::vector<const double *>(size(), &REFERENCE_VALUE));
}

void Household::assign_rate_per_contact(int from, int to, const std::string &parID) { this->_which_rate_[from][to] = parID; }
void Household::set_contact_rate(int from, int to, const double *x) { this->contact_rates[from][to] = x; }
void Household::set_loglike(double x)
{
	this->buffer = this->house_loglike;
	this->house_loglike = x;
}
void Household::restore_loglike() { this->house_loglike = this->buffer; }

//===== GETTER =====//
const std::string &Household::get_ID() const { return this->ID; };
int Household::get_index() const { return this->index; }
int Household::size() const { return this->house_size; }
const std::string &Household::get_group_rate_map(int from, int to) const { return this->_which_rate_[from][to]; }
double Household::get_contact_rate(int from, int to) const { return *this->contact_rates[from][to]; }
bool Household::is_infectious(int i) const { return this->_is_infectious_[i]; }
bool Household::is_susceptible(int i) const { return this->_is_susceptible_[i]; }
bool Household::is_missing_status(int i) const { return this->_is_missing_status_[i]; }
int Household::get_num_infectious() const { return this->num_infectious; }
int Household::get_num_susceptible() const { return this->num_susceptible; }
int Household::get_num_missing_status() const { return this->num_missing_status; }
const std::set<int> &Household::get_infectious_set() const { return this->infectious_set; }
const std::set<int> &Household::get_susceptible_set() const { return this->susceptible_set; }
const std::set<int> &Household::get_missing_status_set() const { return this->missing_status_set; }

double Household::get_loglike() const { return this->house_loglike; }
double Household::get_delta_loglike() const { return this->house_loglike - this->buffer; }

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