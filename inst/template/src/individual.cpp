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

Individual::Individual() {}
Individual::Individual(const std::string &ID_, int status_, const dictargs &data_)
{
	// associate the new ID
	this->ID = ID_;

	// associate the new status
	if (status_ < MISSING_STATUS)
	{
		this->status = MISSING_STATUS;

		std::ostringstream oss;
		oss << "Status '" << status_ << "' of individual " << this->ID << " interpreted as missing."
			<< " Use (=)" << SUSCEPTIBLE << " for susceptible, (>=)" << INFECTIOUS << " for infectious,"
			<< " and (=) " << MISSING_STATUS << " for missing status.";
		tools::warning(oss.str());
	}
	else
	{
		this->status = status_;
	}

	// associate the new data
	this->data = data_;
}
Individual::~Individual() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Individual : implementations  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== SETTER =====//
void Individual::set_status(int s) { this->status = s; }

void Individual::set_rate(const std::string &rateID, const double *x) { this->rate_ptr.emplace(rateID, x); }

/* Disable data and augmented variable assertions. */
void Individual::disable_assert()
{
	_set_augm_method_ = &Individual::_set_augm_default_;
	_get_augm_method_ = &Individual::_get_augm_default_;
	_get_data_method_ = &Individual::_get_data_default_;
	_get_data_pointer_method_ = &Individual::_get_data_pointer_default_;
}

//===== GETTER =====//
std::string Individual::get_ID() const { return this->ID; }
int Individual::get_status() const { return this->status; }
double Individual::get_data(const std::string &dataID) const { return (this->*_get_data_method_)(dataID); }
const double *Individual::get_data_pointer(const std::string &dataID) const { return (this->*_get_data_pointer_method_)(dataID); }
double Individual::get_group_rate(const std::string &rateID) const { return *this->rate_ptr.at(rateID); }
std::unordered_map<std::string, std::string> Individual::get_group_rate_map() const { return this->_which_rate_; }

//===== AUGMENTED VARIABLES =====//

//---- Getter -----//
double Individual::get_augm(const std::string &augmID) const { return (this->*_get_augm_method_)(augmID); }
double Individual::augm_model_logpdf(const std::string &augmID) const { return this->augm_var.at(augmID).distr_logpdf(); }
double Individual::augm_model_logratio(const std::string &augmID) const { return this->augm_var.at(augmID).distr_logratio(); }
double Individual::augm_propose_logratio(const std::string &augmID) const { return this->augm_var.at(augmID).propose_logratio(); }
double Individual::augm_move_logratio(const std::string &augmID) const { return this->augm_var.at(augmID).get_move_logratio(); }

//---- Setter -----//
void Individual::assign_augm(const std::string &augmID, double value,
							 const std::string &model_name, dictargs_ptr model_args_ptr, dictargs range,
							 const std::string &rule_name, dictargs rule_args, double rule_from)
{
	Base_Variable new_augm(augmID, value);
	new_augm.assign_distribution(model_name, model_args_ptr, range);
	new_augm.assign_proposer(rule_name, rule_args, rule_from);
	this->augm_var[augmID] = new_augm;
}
void Individual::set_augm_value(const std::string &augmID, double x) { this->augm_var.at(augmID).set_value(x); }
void Individual::augm_model_rand(const std::string &augmID, sys_rnd_gen &prng) { this->augm_var.at(augmID).rand(prng); }

/* Propose augmented and return true if valid move, otherwise restore buffer.  */
bool Individual::propose_valid_augm_move(const std::string &augmID, sys_rnd_gen &prng)
{
	Base_Variable &this_augm = augm_var.at(augmID);

	// propose new augmented value and store old
	this_augm.propose(prng); 

	// proceed only if valid
	if (this_augm.valid_move())
	{
		return true;
	}
	else
	{
		// reassign buffered value
		this_augm.restore_value();
		return false;
	}
}

void Individual::augm_restore_value(const std::string &augmID) { this->augm_var.at(augmID).restore_value(); }
void Individual::clear_augm(const std::string &augmID)
{
	Base_Variable new_augm(augmID);
	this->augm_var[augmID] = new_augm;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Individual : data  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the individal has the data defined */
void Individual::assert_data(const std::string &dataID) const
{
	if (data.find(dataID) == data.end())
	{
		throw std::invalid_argument(
			"Unknown data '" + dataID + "' for individual " + ID + ".");
	}
}

/* Default method to get the value of data. */
double Individual::_get_data_default_(const std::string &dataID) const
{
	return data.at(dataID);
}

/* Safe method to get the value of data with previous assertion. */
double Individual::_get_data_assert_(const std::string &dataID) const
{
	assert_data(dataID);
	return _get_data_default_(dataID);
}

/* Safe method to get a data pointer with previous assertion. */
const double *Individual::_get_data_pointer_default_(const std::string &dataID) const
{
	return &(data.find(dataID)->second);
}

/* Safe method to get a data pointer with previous assertion. */
const double *Individual::_get_data_pointer_assert_(const std::string &dataID) const
{
	assert_data(dataID);
	return _get_data_pointer_default_(dataID);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Individual : group rates  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign a model rate r */
void Individual::assign_group_rate(const std::string &rateID, const std::string &parID)
{
	if (_which_rate_.find(rateID) == _which_rate_.end())
	{
		// if rateID was undeclared, add it
		_which_rate_.emplace(rateID, parID);
	}
	else
	{
		throw std::invalid_argument(
			"Multiple definitions of group rate '" + rateID + "' for individual " + ID + ".");
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Individual : augmented variables  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise raise error. */
void Individual::assert_augm(const std::string &augmID) const
{
	if (augm_var.find(augmID) == augm_var.end())
	{
		throw std::invalid_argument(
			"Unknown augmented variable '" + augmID + "' for individual " + ID + ".");
	}
}

/* Default method to get the value of an augmented variable. */
double Individual::_get_augm_default_(const std::string &augmID) const
{
	return augm_var.at(augmID).get_value();
}

/* Safe method to get the value of an augmented variable with previous assertion. */
double Individual::_get_augm_assert_(const std::string &augmID) const
{
	assert_augm(augmID);
	return _get_augm_default_(augmID);
}

/* Default method to set the value of an augmented variable. */
void Individual::_set_augm_default_(const std::string &augmID, double new_value)
{
	augm_var.at(augmID).set_value(new_value);
}

/* Safe method to set the value of an augmented variable with previous assertion. */
void Individual::_set_augm_assert_(const std::string &augmID, double new_value)
{
	assert_augm(augmID);
	_set_augm_default_(augmID, new_value);
}