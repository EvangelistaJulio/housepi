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
///////////////////////////////////  Individual : constructors  ///////////////////////////////////
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
void Individual::set_status(int s)
{
	this->status = s;
}

void Individual::set_rate(const std::string &rateID, const double *x)
{
	this->rate_ptr.emplace(rateID, x);
}

/* Disable data and augmented variable assertions. */
void Individual::disable_assert()
{
	this->_func_set_augm_ = &Individual::_set_augm_default_;
	this->_func_get_augm_ = &Individual::_get_augm_default_;
	this->_func_get_data_ = &Individual::_get_data_default_;
	this->_func_get_data_pointer_ = &Individual::_get_data_pointer_default_;
}

//===== GETTER =====//
const std::string &Individual::get_ID() const { return this->ID; }
int Individual::get_status() const { return this->status; }
double Individual::get_data(const std::string &dataID) const { return (this->*_func_get_data_)(dataID); }
const double *Individual::get_data_pointer(const std::string &dataID) const { return (this->*_func_get_data_pointer_)(dataID); }
double Individual::get_group_rate(const std::string &rateID) const { return *this->rate_ptr.at(rateID); }
const std::unordered_map<std::string, std::string> &Individual::get_group_rate_map() const { return this->_which_rate_; }

//===== AUGMENTED VARIABLES =====//

//---- Getter -----//
double Individual::get_augm(const std::string &augmID) const { return (this->*_func_get_augm_)(augmID); }
double Individual::augm_model_logpdf(const std::string &augmID) const { return this->augm_map.at(augmID).distr_logpdf(); }
double Individual::augm_model_logratio(const std::string &augmID) const { return this->augm_map.at(augmID).distr_logratio(); }
double Individual::augm_propose_logratio(const std::string &augmID) const { return this->augm_map.at(augmID).propose_logratio(); }
double Individual::augm_move_logratio(const std::string &augmID) const { return this->augm_map.at(augmID).get_move_logratio(); }

//---- Setter -----//
void Individual::assign_augm(const std::string &augmID, double value, bool is_discr,
							 const std::string &model_name, dictargs_ptr model_args_ptr, dictargs range,
							 const std::string &rule_name, double rule_scale, double rule_from)
{
	this->augm_map[augmID] = Base_Variable(augmID, value, is_discr);
	this->augm_map[augmID].assign_distribution(model_name, model_args_ptr, range);
	this->augm_map[augmID].assign_proposer(rule_name, rule_scale, rule_from);
}

/* Clear the augmented variable definitions. */
void Individual::clear_augm(const std::string &augmID)
{
	this->augm_map[augmID] = Base_Variable(augmID);
}

/* Set a value for the augmented variable. */
void Individual::set_augm_value(const std::string &augmID, double x)
{
	this->augm_map.at(augmID).set_value(x);
}

/* Draw a value from the augmented variable distribution. */
void Individual::augm_model_rand(const std::string &augmID, sys_rnd_gen &prng)
{
	this->augm_map.at(augmID).rand(prng);
}

/* Propose augmented and return true if valid move, otherwise restore buffer. */
bool Individual::propose_valid_augm_move(const std::string &augmID, sys_rnd_gen &prng)
{
	Base_Variable &this_augm = augm_map.at(augmID);

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

/* Restore buffer value of augmented variable. */
void Individual::restore_augm_value(const std::string &augmID)
{
	this->augm_map.at(augmID).restore_value();
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
	if (this->data.find(dataID) == this->data.end())
	{
		throw std::invalid_argument(
			"Unknown data '" + dataID + "' for individual " + ID + ".");
	}
}

/* Default method to get the value of data. */
double Individual::_get_data_default_(const std::string &dataID) const
{
	return this->data.at(dataID);
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
	return &(this->data.find(dataID)->second);
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
	if (this->_which_rate_.find(rateID) == this->_which_rate_.end())
	{
		// if rateID was undeclared, add it
		this->_which_rate_.emplace(rateID, parID);
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
	if (this->augm_map.find(augmID) == this->augm_map.end())
	{
		throw std::invalid_argument(
			"Unknown augmented variable '" + augmID + "' for individual " + ID + ".");
	}
}

/* Default method to get the value of an augmented variable. */
double Individual::_get_augm_default_(const std::string &augmID) const
{
	return this->augm_map.at(augmID).get_value();
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
	this->augm_map.at(augmID).set_value(new_value);
}

/* Safe method to set the value of an augmented variable with previous assertion. */
void Individual::_set_augm_assert_(const std::string &augmID, double new_value)
{
	assert_augm(augmID);
	_set_augm_default_(augmID, new_value);
}
