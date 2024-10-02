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
	// associate index
	if (new_data.find("index") == new_data.end())
	{
		basics::raise_error("Missing index for Individual object.");
	}
	else
	{
		set_index(new_data.at("index"));
	}

	// associate status if provided
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
			oss << "Status '" << new_status << "' of individual '" << get_index() << "' interpreted as missing."
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
		oss << "Unknown data '" << dataID << "' for individual " << get_index() << ".";
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
void Individual::assign_rate_per_group(const string &rateID, const string &parID)
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
			<< "' for individual " << get_index() << "."
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
			<< "' for individual " << get_index() << "."
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
	augm_var.at(augmID).valid(new_value);
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

/* Assign an augmented variable to the individual. */
void Individual::assign_augm(const string &augmID, double value,
							 const string &model_name, const dictargs_ptr &model_args_ext,
							 const string &rule_name, const dictargs &rule_args, const dictargs &range, double rule_from)
{
	// check that lists have equal sizes
	if (model_args_ext.find(IO_MISMATCH) != model_args_ext.end())
	{
		raise_error("Model arguments, values and types of augmented variable '" + augmID + "' should have equal size.");
	}
	if (rule_args.find(IO_MISMATCH) != rule_args.end())
	{
		raise_error("Rule arguments, values and types of augmented variable '" + augmID + "' should have equal size.");
	}

	Augmented new_augm(value, model_name, model_args_ext,
					   rule_name, rule_args, range, rule_from);

	// overwrite augmented variable if it has multiple definitions
	// FIXME : is this relevant?
	if (augm_var.find(augmID) == augm_var.end())
	{
		augm_var.emplace(augmID, new_augm);
	}
	else
	{
		augm_var.at(augmID) = new_augm;
	}
}

/* Init augmented variable within a random initial value */
void Individual::init_augm(const string &augmID, sys_rnd_gen &prng)
{
	double init_value = augm_model_rand(augmID, prng);
	int num_attempts = 0;

	while ((!augm_valid(augmID, init_value)))
	{
		if (num_attempts > MAX_NUMBER_ATTEMPTS)
		{
			std::ostringstream oss;
			oss << "Exceeded attempts  (" << ">" << MAX_NUMBER_ATTEMPTS << ")"
				<< " in random initialization of '" << augmID << "'"
				<< " for individual " << get_index() << ".";
			raise_error(oss.str());
		}
		init_value = augm_model_rand(augmID, prng);
		++num_attempts;
	}
	set_augm(augmID, init_value);
}

/* Propose a new value for the augmented variable. */
double Individual::augm_propose(const string &augmID, sys_rnd_gen &prng)
{
	return augm_var.at(augmID).propose(prng);
}

/* Log ratio of the proposal rule associated to the move old_value to new_value. */
double Individual::augm_propose_logratio(const string &augmID, double old_value, double new_value) const
{
	return augm_var.at(augmID).propose_logratio(old_value, new_value);
}

/* Log ratio of the model associated to the move old_value to new_value.  */
double Individual::augm_model_logratio(const string &augmID, double old_value, double new_value) const
{
	return augm_var.at(augmID).model_logratio(old_value, new_value);
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
		oss << "Unknown data '" << dataID << "' for household " << get_index() << ".";
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
		oss << "Impossibile to assign individual '" << this_indiv.get_index()
			<< "' : unknown status '" << status << "'.";
		raise_error(oss.str());
	}

	//----- add the new individual object -----//
	indiv.push_back(this_indiv);

	//----- update household size -----//
	++house_size;
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
			_is_infectious_[i] = false;
		}
		else if (old_status == SUSCEPTIBLE)
		{
			num_susceptible--;
			susceptible_set.erase(i);
			_is_susceptible_[i] = false;
		}
		else if (old_status == MISSING_STATUS)
		{
			num_missing_status--;
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
	dgn_loglike = vector<double>(sample_size, std::numeric_limits<double>::quiet_NaN());
}

/* update the values for the Bayesian diagnostics */
void Household::update_diagnostics()
{
	dgn_loglike[dgn_count++] = house_loglike;
}