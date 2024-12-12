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

void Individual::assign_all_group_rates(const dict_rates_ptr &group_rate_ptrs)
{
	for (const auto &pair : group_rate_ptrs)
	{
		std::string rateID = pair.first;
		std::vector<const double *> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate '" + rateID + "' for individual " + ID + ".");
		else
			this->rate_ptr.emplace(rateID, vec_par_ptr);
	}
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
double Individual::get_data(const dataID &d) const { return (this->*_func_get_data_)(d); }
const double *Individual::get_data_pointer(const dataID &d) const { return (this->*_func_get_data_pointer_)(d); }
double Individual::get_group_rate(const std::string &rateID) const
{
	return *this->rate_ptr.at(rateID)[0]; // FIXME : <== DEBUGGING HERE
}

//===== AUGMENTED VARIABLES =====//

//---- Getter -----//
double Individual::get_augm(const augmID &a) const { return (this->*_func_get_augm_)(a); }
double Individual::augm_model_logpdf(const augmID &a) const { return this->augm_map.at(a).distr_logpdf(); }
double Individual::augm_model_logratio(const augmID &a) const { return this->augm_map.at(a).distr_logratio(); }
double Individual::augm_propose_logratio(const augmID &a) const { return this->augm_map.at(a).propose_logratio(); }
double Individual::augm_move_logratio(const augmID &a) const { return this->augm_map.at(a).get_move_logratio(); }

//---- Setter -----//
void Individual::assign_augm(const augmID &a, double value, bool is_discr,
							 const std::string &model_name, dictargs_ptr model_args_ptr, dictargs range,
							 const std::string &rule_name, double rule_scale, double rule_from)
{
	this->augm_map[a] = Base_Variable(a, value, is_discr);
	this->augm_map[a].assign_distribution(model_name, model_args_ptr, range);
	this->augm_map[a].assign_proposer(rule_name, rule_scale, rule_from);
}

/* Clear the augmented variable definitions. */
void Individual::clear_augm(const augmID &a)
{
	this->augm_map[a] = Base_Variable(a);
}

/* Set a value for the augmented variable. */
void Individual::set_augm_value(const augmID &a, double x)
{
	this->augm_map.at(a).set_value(x);
}

/* Draw a value from the augmented variable distribution. */
void Individual::augm_model_rand(const augmID &a, sys_rnd_gen &prng)
{
	this->augm_map.at(a).rand(prng);
}

/* Propose augmented and return true if valid move, otherwise restore buffer. */
bool Individual::propose_valid_augm_move(const augmID &a, sys_rnd_gen &prng)
{
	Base_Variable &this_augm = augm_map.at(a);

	this_augm.propose(prng); // propose new value and buffer old

	if (this_augm.valid_move()) // proceed only if valid
		return true;
	else
		this_augm.restore_value(); // reassign buffered value

	return false;
}

/* Restore buffer value of augmented variable. */
void Individual::restore_augm_value(const augmID &a)
{
	this->augm_map.at(a).restore_value();
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
bool Individual::assert_data(const dataID &d) const
{
	if (this->data.find(d) == this->data.end())
		throw std::invalid_argument(
			"Unknown data '" + d + "' for individual " + ID + ".");
	else
		return true;

	return false;
}

/* Default method to get the value of data. */
double Individual::_get_data_default_(const dataID &d) const
{
	return this->data.at(d);
}

/* Safe method to get the value of data with previous assertion. */
double Individual::_get_data_assert_(const dataID &d) const
{
	if (assert_data(d))
		return _get_data_default_(d);
	else
		return INF_DOUBLE;
}

/* Safe method to get a data pointer with previous assertion. */
const double *Individual::_get_data_pointer_default_(const dataID &d) const
{
	return &(this->data.find(d)->second);
}

/* Safe method to get a data pointer with previous assertion. */
const double *Individual::_get_data_pointer_assert_(const dataID &d) const
{
	if (assert_data(d))
		return _get_data_pointer_default_(d);
	else
		return nullptr;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Individual : augmented variables  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the given parameter is defined, otherwise raise error. */
bool Individual::assert_augm(const augmID &a) const
{
	if (this->augm_map.find(a) == this->augm_map.end())
		throw std::invalid_argument(
			"Unknown augmented variable '" + a + "' for individual " + ID + ".");
	else
		return true;

	return false;
}

/* Default method to get the value of an augmented variable. */
double Individual::_get_augm_default_(const augmID &a) const
{
	return this->augm_map.at(a).get_value();
}

/* Safe method to get the value of an augmented variable with previous assertion. */
double Individual::_get_augm_assert_(const augmID &a) const
{
	if (assert_augm(a))
		return _get_augm_default_(a);
	else
		return INF_DOUBLE;
}

/* Default method to set the value of an augmented variable. */
void Individual::_set_augm_default_(const augmID &a, double new_value)
{
	this->augm_map.at(a).set_value(new_value);
}

/* Safe method to set the value of an augmented variable with previous assertion. */
void Individual::_set_augm_assert_(const augmID &a, double new_value)
{
	if (assert_augm(a))
		_set_augm_default_(a, new_value);
}
