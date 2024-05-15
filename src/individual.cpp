// individual.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "individual.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////// Individual : constructor ////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Constructor with a database row. */
Individual::Individual(sdargs new_data)
{
	index = new_data.at("index");
	status = new_data.at("status");
	data = new_data;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Individual : methods //////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Return a feature of the Individual as provided in the database */
double Individual::get_data(std::string dataID) const
{
	return data.at(dataID);
}

/* Check if the individal has the data defined */
void Individual::assert_data(std::string dataID)
{
	if (data.find(dataID) == data.end())
	{
		std::ostringstream oss;
		oss << "Unknown data `" << dataID << "` for individual " << get_index() << ".";
		display_error_message(oss);
		std::exit(EXIT_FAILURE);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////// Individual : methods (Groutp Rates) /////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign a model rate */
void Individual::assign_rate_per_group(std::string rateID, std::string parID)
{
	if (_which_rate.find(rateID) == _which_rate.end())
	{
		// if rateID was undeclared, add it
		_which_rate.emplace(rateID, parID);
	}
	else
	{
		// rateID was previously declared
		std::ostringstream oss;
		oss << "Individual " << get_index() << " has multiple definitions for group rate `" << rateID << "`.";
		display_error_message(oss);
		std::exit(EXIT_FAILURE);
	}
}

/* Link the rate value with a pointer */
void Individual::_set_rate(std::string rateID, const double *x)
{
	_rate_ptr.emplace(rateID, x);
}

/* Return value of model rate */
double Individual::model_rate(std::string rateID) const
{
	return *_rate_ptr.at(rateID);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////// Individual : methods (Augmented Variables) /////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Init augmented variable within a random initial value */
void Individual::init_augm(std::string augmID, std::mt19937_64 &rng)
{
	int num_attempts = 0;
	double init_value = augm_model_rand(augmID, rng);

	while ((!augm_valid(augmID, init_value)))
	{
		if (num_attempts >= MAX_NUMBER_ATTEMPTS)
		{
			std::ostringstream oss;
			oss << "Number of attempts exceeded (" << ">=" << MAX_NUMBER_ATTEMPTS << ")"
				<< " in random initialization of `" << augmID << "`"
				<< " for individual " << get_index();
			display_error_message(oss);
			std::exit(EXIT_FAILURE);
		}
		init_value = augm_model_rand(augmID, rng);
		num_attempts++;
	}
	set_augm(augmID, init_value);
}

/* Add an augmented variable to the individual. */
void Individual::add_augm(std::string augmID, Augmented new_augm)
{
	if (augm_var.find(augmID) == augm_var.end())
	{
		augm_var.emplace(augmID, new_augm);
	}
	else
	{
		std::ostringstream oss;
		oss << "Individual " << get_index() << " has multiple definitions for augmented variable `" << augmID << "`.";
		display_error_message(oss);
		std::exit(EXIT_FAILURE);
	}
}

/* Propose a new value for the augmented variable. */
double Individual::augm_propose(std::string augmID, std::mt19937_64 &prng)
{
	return augm_var.at(augmID).propose(prng);
}

/* Log ratio of the proposal rule associated to the move old_value to new_value. */
double Individual::augm_propose_logratio(std::string augmID, double old_value, double new_value)
{
	return augm_var.at(augmID).propose_logratio(old_value, new_value);
}

/* Log ratio of the model associated to the move old_value to new_value.  */
double Individual::augm_model_logratio(std::string augmID, double old_value, double new_value)
{
	return augm_var.at(augmID).model_logratio(old_value, new_value);
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////// Household : methods ////////////////////////////////////// //
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
		display_error_message(oss);
		std::exit(EXIT_FAILURE);
	}
}

/* Return a feature of the Household which must be provided during initialization. */
double Household::get_data(std::string dataID) const
{
	return house_data.at(dataID);
}

/* Link an individual to the household */
void Household::add_individual(Individual this_indiv)
{
	int status = this_indiv.get_status();
	indiv.push_back(this_indiv);

	// status list
	if (status == INFECTIOUS)
	{
		infectious_list.push_back(house_size);
	}
	else if (status == SUSCEPTIBLE)
	{
		susceptibles_list.push_back(house_size);
	}
	else
	{
		unknown_status_list.push_back(house_size);
	}
	house_size++;
}

/* Return the indexes of the infectious, then the susceptibles, unknown statuses */
// FIXME : is this necessary anymore?
std::vector<int> Household::get_house_members()
{
	std::vector<int> members;
	members.reserve(house_size); // Reserve space

	// Copy infectious
	members.insert(members.end(), infectious_list.begin(), infectious_list.end());

	// Copy susceptibles
	members.insert(members.end(), susceptibles_list.begin(), susceptibles_list.end());

	// Copy unknown
	members.insert(members.end(), unknown_status_list.begin(), unknown_status_list.end());

	return members;
}

/* Change_status of individual i to new status (without controls).
	WARNING : this introduces augmented variables -> it requires RJMCMC in case of inference.*/
void Household::change_status(int i, int new_status)
{
	int old_status = indiv[i].get_status();

	if (old_status != new_status)
	{
		// erase individual i from previous list
		if (old_status == SUSCEPTIBLE)
		{
			susceptibles_list.erase(find(susceptibles_list.begin(), susceptibles_list.end(), i));
		}
		else if (old_status == INFECTIOUS)
		{
			infectious_list.erase(find(infectious_list.begin(), infectious_list.end(), i));
		}
		else
		{
			unknown_status_list.erase(find(unknown_status_list.begin(), unknown_status_list.end(), i));
		}

		// add individual i to new list
		if (new_status == SUSCEPTIBLE)
		{
			susceptibles_list.push_back(i);
		}
		else if (new_status == INFECTIOUS)
		{
			infectious_list.push_back(i);
		}
		else
		{
			unknown_status_list.push_back(i);
		}

		// change the status of individual to new
		indiv[i].set_status(new_status);
	}
}

/* Set to zero the counts of given augmented variables moves. */
void Household::reset_augm_accept_rate(std::string augmID)
{
	num_accept_augm_moves[augmID] = 0;
	num_propos_augm_moves[augmID] = 0;
}

/* Update counts for acceptance rate of given augemented variable. */
void Household::update_augm_accept_rate(std::string augmID, bool is_accepted)
{
	num_accept_augm_moves[augmID] += is_accepted;
	num_propos_augm_moves[augmID]++;
}

/* update the values for the Bayesian diagnostics */
void Household::update_diagnostics()
{
	dgn_sum_likelihood += std::exp(loglike);
	dgn_sum_loglike += loglike;
	dgn_sum_sqr_loglike += std::pow(loglike, 2);
	dgn_ss++;
}
