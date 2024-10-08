/*
 * individual.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "individual.h"

using basics::raise_error;

//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

bool Unknown_Status_Warning = false; // unknown status

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
Individual::Individual(dictargs new_data)
{
	// check that an index has been associated by call to build
	if (new_data.find("index") == new_data.end())
	{
		basics::raise_error("Missing index for Individual object.");
	}

	// associate the status if provided
	if (new_data.find("status") == new_data.end())
	{
		set_status(MISSING_STATUS);
	}
	else
	{
		int new_status = new_data.at("status");

		//----- Check for the validity of the status -----//
		if (VALID_STATUS.find(new_status) == VALID_STATUS.end())
		{
			set_status(MISSING_STATUS);

			std::ostringstream oss;
			oss << "Status '" << new_status << "' of individual " << get_data("index") << " interpreted as missing."
				<< " Use " << SUSCEPTIBLE << "/" << INFECTIOUS << " for susceptible/infectious status."
				<< " Use " << MISSING_STATUS << " for missing status.";
			basics::raise_warning(oss.str(), Unknown_Status_Warning);
		}
		else
		{
			status = new_status;
		}
	}

	data = new_data;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Individual : methods (general)  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the individal has the data defined */
void Individual::assert_data(const string &dataID) const
{
	if (data.find(dataID) == data.end())
	{
		std::ostringstream oss;
		oss << "Unknown data '" << dataID << "' for individual " << get_data("index") << ".";
		oss << " Check your model and your data definitions.";
		raise_error(oss.str());
	}
}

/* Default method to get the value of data. */
double Individual::_get_data_default_(const string &dataID) const
{
	return data.at(dataID);
}

/* Safe method to get the value of data with previous assertion. */
double Individual::_get_data_assert_(const string &dataID) const
{
	assert_data(dataID);
	return _get_data_default_(dataID);
}

/* Safe method to get a data pointer with previous assertion. */
const double *Individual::_get_data_pointer_default_(const string &dataID) const
{
	return &(data.find(dataID)->second);
}

/* Safe method to get a data pointer with previous assertion. */
const double *Individual::_get_data_pointer_assert_(const string &dataID) const
{
	assert_data(dataID);
	return _get_data_pointer_default_(dataID);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Individual : methods (group rates)  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign a model rate r */
void Individual::assign_group_rate(const string &rateID, const string &parID)
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
		oss << "Multiple definitions of group rate '" << rateID
			<< "' for individual " << get_data("index") << "."
			<< " Check your model.";
		raise_error(oss.str());
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Individual : methods (augmented variables)  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise raise error. */
void Individual::assert_augm(const string &augmID) const
{
	if (augm_var.find(augmID) == augm_var.end())
	{
		std::ostringstream oss;
		oss << "Unknown augmented variable '" << augmID
			<< "' for individual " << get_data("index") << "."
			<< " Check your model and your augmented variable definitions.";
		raise_error(oss.str());
	}
}

/* Default method to get the value of an augmented variable. */
double Individual::_get_augm_default_(const string &augmID) const
{
	return augm_var.at(augmID).get_value();
}

/* Safe method to get the value of an augmented variable with previous assertion. */
double Individual::_get_augm_assert_(const string &augmID) const
{
	assert_augm(augmID);
	return _get_augm_default_(augmID);
}

/* Default method to set the value of an augmented variable. */
void Individual::_set_augm_default_(const string &augmID, double new_value)
{
	augm_var.at(augmID).set_value(new_value);
}

/* Safe method to set the value of an augmented variable with previous assertion. */
void Individual::_set_augm_assert_(const string &augmID, double new_value)
{
	assert_augm(augmID);
	_set_augm_default_(augmID, new_value);
}

/* Disable data and augmented variable assertions. */
void Individual::disable_assert()
{
	_set_augm_method_ = &Individual::_set_augm_default_;
	_get_augm_method_ = &Individual::_get_augm_default_;
	_get_data_method_ = &Individual::_get_data_default_;
	_get_data_pointer_method_ = &Individual::_get_data_pointer_default_;
}

/* Propose augmented and return true if valid move, otherwise restore buffer.  */
bool Individual::propose_valid_augm_move(const string &augmID, sys_rnd_gen &prng)
{
	// propose new augmented value and store old in buffer
	augm_propose(augmID, prng);

	// proceed only if valid
	if (augm_valid_move(augmID))
	{
		return true;
	}
	else
	{
		// reassign buffered value
		augm_restore_value(augmID);
		return false;
	}
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
void Household::assign_data(const string &dataID, double value)
{
	if (house_data.find(dataID) == house_data.end())
	{
		// if rateID was undeclared, add it
		house_data.emplace(dataID, value);
	}
	else
	{
		// rateID was previously declared, overwrite it
		house_data.at(dataID) = value;
	}
}

/* Return a feature of the Household which must be provided during initialization. */
double Household::get_data(const string &dataID) const
{
	assert_data(dataID);
	return house_data.at(dataID);
}

/* Check if the data is defined for the Household. */
void Household::assert_data(const string &dataID) const
{
	if (house_data.find(dataID) == house_data.end())
	{
		std::ostringstream oss;
		oss << "Unknown data '" << dataID << "' for household " << get_data("index") << ".";
		raise_error(oss.str());
	}
}

/* Link an individual to the household */
void Household::assign_individual(Individual this_indiv)
{
	int status = this_indiv.get_status();

	//----- update status lists within the household -----//
	if (status == INFECTIOUS)
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
		std::ostringstream oss;
		oss << "Impossibile to assign individual '" << this_indiv.get_data("index")
			<< "' : unknown status '" << status << "'.";
		raise_error(oss.str());
	}

	//----- add the new individual object -----//
	indiv.push_back(this_indiv);

	//----- update household size -----//
	++house_size;
}

/* Change_status of individual i to new status (without controls).
	NOTE : this introduces augmented variables -> it requires RJMCMC in case of inference. */
bool Household::change_status(int i, int new_status)
{
	int old_status = indiv[i].get_status();

	if (new_status == old_status)
	{
		return false;
	}
	else
	{
		if (old_status == INFECTIOUS)
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
			oss << "Impossibile to change status: unknown old status '" << old_status << "'.";
			raise_error(oss.str());
		}

		if (new_status == INFECTIOUS)
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
			oss << "Impossibile to change status: unknown new status '" << new_status << "'.";
			raise_error(oss.str());
		}
		indiv[i].set_status(new_status);

		return true;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////  Household : methods (monte carlo)  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Set to zero the counts of given augmented variables moves. */
void Household::reset_augm_accept_rate(const string &augmID)
{
	num_accept_augm_moves[augmID] = 0;
	num_propos_augm_moves[augmID] = 0;
}

/* Update counts for acceptance rate of given augemented variable. */
void Household::update_augm_accept_rate(const string &augmID, bool is_accepted)
{
	num_accept_augm_moves[augmID] += is_accepted;
	num_propos_augm_moves[augmID]++;
}

/* update the values for the Bayesian diagnostics */
void Household::init_diagnostics(ulong sample_size)
{
	dgn_count = 0;
	dgn_loglike.clear();
	dgn_loglike = vector<double>(sample_size);
}

/* update the values for the Bayesian diagnostics */
void Household::update_diagnostics()
{
	dgn_loglike[dgn_count++] = house_loglike;
}