/*
 * individual.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "individual.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Individual : constructor  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Constructor with a database row. */
Individual::Individual(sdargs new_data)
{
	set_index(new_data.at("index"));

	// associate status if provided
	if (new_data.find("status") == new_data.end())
	{
		status = MISSING_STATUS;
	}
	else
	{
		status = new_data.at("status");
	}

	data = new_data;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Individual : methods (general)  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the individal has the data defined */
void Individual::assert_data(std::string d) const
{
	if (data.find(d) == data.end())
	{
		std::ostringstream oss;
		oss << "Unknown data `" << d << "` for individual " << get_index() << ".";
		screen::raise_error(oss.str());
	}
}

/* Default method to get the value of data. */
double Individual::_get_data_default_(const std::string d) const
{
	return data.at(d);
}

/* Safe method to get the value of data with previous assertion. */
double Individual::_get_data_assert_(const std::string d) const
{
	assert_data(d);
	return _get_data_default_(d);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Individual : methods (group rates)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign a model rate */
void Individual::assign_rate_per_group(std::string rateID, std::string parID)
{
	if (which_rate.find(rateID) == which_rate.end())
	{
		// if rateID was undeclared, add it
		which_rate.emplace(rateID, parID);
	}
	else
	{
		// rateID was previously declared
		std::ostringstream oss;
		oss << "Individual " << get_index() << " has multiple definitions for group rate `" << rateID << "`.";
		screen::raise_error(oss.str());
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Individual : methods (augmented variables)  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise raise error. */
void Individual::assert_augm(const std::string a) const
{
	if (augm_var.find(a) == augm_var.end())
	{
		std::ostringstream oss;
		oss << "Unknown augmented variable `" << a << "` for individual " << get_index() << ".";
		screen::raise_error(oss.str());
	}
}

/* Default method to get the value of an augmented variable. */
double Individual::_get_augm_default_(const std::string a) const
{
	return augm_var.at(a).get_value();
}

/* Safe method to get the value of an augmented variable with previous assertion. */
double Individual::_get_augm_assert_(const std::string a) const
{
	assert_augm(a);
	return _get_augm_default_(a);
}

/* Disable data and augmented variable assertions. */
void Individual::disable_assert()
{
	_get_augm_method_ = &Individual::_get_augm_default_;
	_get_data_method_ = &Individual::_get_data_default_;
}

/* Assign an augmented variable to the individual. */
void Individual::assign_augm(const std::string a, Augmented new_augm)
{
	if (augm_var.find(a) == augm_var.end())
	{
		augm_var.emplace(a, new_augm);
	}
	else
	{
		std::ostringstream oss;
		oss << "Individual " << get_index() << " has multiple definitions for augmented variable `" << a << "`.";
		screen::raise_error(oss.str());
	}
}

/* Init augmented variable within a random initial value */
void Individual::init_augm(const std::string a, sys_rnd_gen &rng)
{
	int num_attempts = 0;
	double init_value = augm_model_rand(a, rng);

	while ((!augm_valid(a, init_value)))
	{
		if (num_attempts >= MAX_NUMBER_ATTEMPTS)
		{
			std::ostringstream oss;
			oss << "Number of attempts exceeded (" << ">=" << MAX_NUMBER_ATTEMPTS << ")"
				<< " in random initialization of `" << a << "`"
				<< " for individual " << get_index();
			screen::raise_error(oss.str());
		}
		init_value = augm_model_rand(a, rng);
		num_attempts++;
	}
	set_augm(a, init_value);
}

/* Propose a new value for the augmented variable. */
double Individual::augm_propose(const std::string a, sys_rnd_gen &prng)
{
	return augm_var.at(a).propose(prng);
}

/* Log ratio of the proposal rule associated to the move old_value to new_value. */
double Individual::augm_propose_logratio(const std::string a, double old_value, double new_value)
{
	return augm_var.at(a).propose_logratio(old_value, new_value);
}

/* Log ratio of the model associated to the move old_value to new_value.  */
double Individual::augm_model_logratio(const std::string a, double old_value, double new_value)
{
	return augm_var.at(a).model_logratio(old_value, new_value);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Household : methods (general)  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign a data feature to the household */
void Household::assign_data(std::string dataID, double value)
{
	if (house_data.find(dataID) == house_data.end())
	{
		// if rateID was undeclared, add it
		house_data.emplace(dataID, value);
	}
	else
	{
		// rateID was previously declared
		std::ostringstream oss;
		oss << "Household " << get_index() << " has multiple definitions for data `" << dataID << "`.";
		screen::raise_error(oss.str());
	}
}

/* Return a feature of the Household which must be provided during initialization. */
double Household::get_data(std::string dataID) const
{
	assert_data(dataID);
	return house_data.at(dataID);
}

/* Check if the data is defined for the Household. */
void Household::assert_data(std::string dataID) const
{
	if (house_data.find(dataID) == house_data.end())
	{
		std::ostringstream oss;
		oss << "Unknown data `" << dataID << "` for household " << get_index() << ".";
		screen::raise_error(oss.str());
	}
}

/* Link an individual to the household */
void Household::assign_individual(Individual this_indiv)
{
	int status = this_indiv.get_status();

	//----- update status lists within the household -----//
	if (status == INFECTIOUS)
	{
		num_infectious++;
		infectious_set.insert(house_size);

		is_infectious.push_back(true);
		is_susceptible.push_back(false);
		is_missing_status.push_back(false);
	}
	else if (status == SUSCEPTIBLE)
	{
		num_susceptible++;
		susceptible_set.insert(house_size);

		is_infectious.push_back(false);
		is_susceptible.push_back(true);
		is_missing_status.push_back(false);
	}
	else if (status == MISSING_STATUS)
	{
		num_missing_status++;
		missing_status_set.insert(house_size);

		is_infectious.push_back(false);
		is_susceptible.push_back(false);
		is_missing_status.push_back(true);
	}
	else
	{
		std::ostringstream oss;
		oss << "Impossibile to assign individual " << this_indiv.get_index() << ": unknown status " << status << "`.";
		screen::raise_error(oss.str());
	}

	//----- add the new individual object -----//
	indiv.push_back(this_indiv);

	//----- update household size -----//
	house_size++;
}

/* Change_status of individual i to new status (without controls).
	WARNING : this introduces augmented variables -> it requires RJMCMC in case of inference.*/
void Household::change_status(int i, int new_status)
{
	int old_status = indiv[i].get_status();

	if (new_status != old_status)
	{
		if (old_status == INFECTIOUS)
		{
			num_infectious--;
			infectious_set.erase(i);
			is_infectious[i] = false;
		}
		else if (old_status == SUSCEPTIBLE)
		{
			num_susceptible--;
			susceptible_set.erase(i);
			is_susceptible[i] = false;
		}
		else if (old_status == MISSING_STATUS)
		{
			num_missing_status--;
			missing_status_set.erase(i);
			is_missing_status[i] = false;
		}
		else
		{
			std::ostringstream oss;
			oss << "Impossibile to change status: unknown old status `" << old_status << "`.";
			screen::raise_error(oss.str());
		}

		if (new_status == INFECTIOUS)
		{
			num_infectious++;
			infectious_set.insert(i);
			is_infectious[i] = true;
		}
		else if (new_status == SUSCEPTIBLE)
		{
			num_susceptible++;
			susceptible_set.insert(i);
			is_susceptible[i] = true;
		}
		else if (new_status == MISSING_STATUS)
		{
			num_missing_status++;
			missing_status_set.insert(i);
			is_missing_status[i] = true;
		}
		else
		{
			std::ostringstream oss;
			oss << "Impossibile to change status: unknown new status `" << new_status << "`.";
			screen::raise_error(oss.str());
		}

		indiv[i].set_status(new_status);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Household : methods (monte carlo)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set to zero the counts of given augmented variables moves. */
void Household::reset_augm_accept_rate(std::string a)
{
	num_accept_augm_moves[a] = 0;
	num_propos_augm_moves[a] = 0;
}

/* Update counts for acceptance rate of given augemented variable. */
void Household::update_augm_accept_rate(std::string a, bool is_accepted)
{
	num_accept_augm_moves[a] += is_accepted;
	num_propos_augm_moves[a]++;
}

/* update the values for the Bayesian diagnostics */
void Household::init_diagnostics(int sample_size)
{
	dgn_ss = sample_size;
	dgn_like = 0.;
	dgn_avg = 0.;
	dgn_sqr = 0.;
}

/* update the values for the Bayesian diagnostics */
void Household::update_diagnostics()
{
	dgn_like += std::exp(house_loglike) / dgn_ss;
	dgn_avg += house_loglike / dgn_ss;
	dgn_sqr += std::pow(house_loglike, 2) / dgn_ss;
}