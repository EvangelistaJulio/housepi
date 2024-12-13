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
/////////////////////////////////////  Group Rate : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== CONSTRUCTOR	=====//

GroupRate::GroupRate() {};
GroupRate::GroupRate(const std::vector<const double *> &pointers)
{
	if (pointers.empty())
		throw std::invalid_argument(
			"No parameters to associate to group rate.");
	else if (pointers.size() == 1)
		_func_getter_ = &_getter_default_;
	else
		_func_getter_ = &_getter_multiple_;

	this->pointers = pointers;
}
GroupRate::~GroupRate() = default;

//===== GETTER =====//
double GroupRate::get() const { (this->*_func_getter_)(); }

double GroupRate::_getter_error_() const
{
	throw std::invalid_argument(
		"Attempt to call getter of undefined group rate.");
	return INF_DOUBLE;
}

double GroupRate::_getter_default_() const { return *this->pointers[0]; }

double GroupRate::_getter_multiple_() const
{
	double output = 1.0;
	for (const double *ptr : pointers)
		output *= *ptr;
	return output;
}

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
Individual::Individual(const std::string &ID,
					   int status,
					   const dataIDmap<double> &data)
{
	// associate the new ID
	this->ID = ID;

	// associate the new status
	if (status < MISSING_STATUS)
	{
		this->status = MISSING_STATUS;

		std::ostringstream oss;
		oss << "Status '" << status << "' of individual " << ID << " interpreted as missing."
			<< " Use (=)" << SUSCEPTIBLE << " for susceptible, (>=)" << INFECTIOUS << " for infectious,"
			<< " and (=) " << MISSING_STATUS << " for missing status.";
		tools::warning(oss.str());
	}
	else
	{
		set_status(status);
	}

	// associate the new data
	this->data = data;
}
Individual::~Individual() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Individual : implementations  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== SETTER =====//
void Individual::set_status(int s) { this->status = s; }

void Individual::assign_all_group_rates(const rateIDmapvec<const double *> &group_rate_ptrs)
{
	for (const auto &pair : group_rate_ptrs)
	{
		rateID r = pair.first;
		std::vector<const double *> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate '" + r + "' for individual " + ID + ".");
		else
			this->rate_ptr[r] = GroupRate(vec_par_ptr);
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
double Individual::get_group_rate(const rateID &r) const { return this->rate_ptr.at(r).get(); }

//===== AUGMENTED VARIABLES =====//

//---- Getter -----//
double Individual::get_augm(const augmID &a) const { return (this->*_func_get_augm_)(a); }
double Individual::augm_model_logpdf(const augmID &a) const { return this->augm_map.at(a).distr_logpdf(); }
double Individual::augm_model_logratio(const augmID &a) const { return this->augm_map.at(a).distr_logratio(); }
double Individual::augm_propose_logratio(const augmID &a) const { return this->augm_map.at(a).propose_logratio(); }
double Individual::augm_move_logratio(const augmID &a) const { return this->augm_map.at(a).get_move_logratio(); }

//---- Setter -----//
void Individual::assign_augm(const augmID &a, double value, bool is_discr,
							 const std::string &model_name, dict_args_ptr model_args_ptr, dict_args range,
							 const std::string &rule_name, double rule_scale, double rule_from)
{
	this->augm_map[a] = Base_Variable(a, value, is_discr);
	this->augm_map[a].assign_distribution(model_name, model_args_ptr, range);
	this->augm_map[a].assign_proposer(rule_name, rule_scale, rule_from);
}

/* Clear the augmented variable definitions. */
void Individual::clear_augm(const augmID &a) { this->augm_map[a] = Base_Variable(a); }

/* Set a value for the augmented variable. */
void Individual::set_augm_value(const augmID &a, double x) { this->augm_map.at(a).set_value(x); }

/* Restore buffer value of augmented variable. */
void Individual::restore_augm_value(const augmID &a) { this->augm_map.at(a).restore_value(); }

/* Propose augmented and return true if valid move, otherwise restore buffer. */
bool Individual::propose_valid_augm_move(const augmID &a, rnd_gen &prng)
{
	Base_Variable &this_augm = augm_map.at(a);

	this_augm.propose(prng); // propose new value and buffer old

	if (this_augm.valid_move()) // proceed only if valid
		return true;
	else
		this_augm.restore_value(); // reassign buffered value

	return false;
}

/* Draw a value from the augmented variable distribution. */
void Individual::augm_model_rand(const augmID &a, rnd_gen &prng) { this->augm_map.at(a).rand(prng); }

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
void Individual::_set_augm_default_(const augmID &a, double xnew)
{
	this->augm_map.at(a).set_value(xnew);
}

/* Safe method to set the value of an augmented variable with previous assertion. */
void Individual::_set_augm_assert_(const augmID &a, double xnew)
{
	if (assert_augm(a))
		_set_augm_default_(a, xnew);
}
