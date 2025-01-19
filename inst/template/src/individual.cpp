/*
 * individual.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "individual.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Individual::Group Rate : class  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== CONSTRUCTOR	=====//
Individual::GroupRate::GroupRate(const std::string &name,
								 const std::vector<const double *> &pointers)
{
	if (!tools::is_na(name))
	{
		this->name = name;

		if (pointers.empty())
			throw std::invalid_argument(
				"No parameters associated to group rate \"" + name + "\".");

		this->n_rates = pointers.size();
		this->pointers = new const double *[pointers.size()];
		for (ushort c = 0; c < this->n_rates; ++c)
			this->pointers[c] = pointers.at(c);
	}
}
Individual::GroupRate::~GroupRate() = default;

//===== GETTER =====//
double Individual::GroupRate::get() const
{
	double output = 1.0;
	for (ushort r = 0; r < n_rates; ++r)
		output *= *pointers[r];

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

// constructor
Individual::Individual(const std::string &ID,
					   int status,
					   const dataIDmap<double> &data)
{
	// associate the new ID
	this->ID = ID;
	// associate the new status
	this->status = (status < MISSING_STATUS) ? MISSING_STATUS : status;
	// associate the new data
	this->data = data;
}

// move constructor // UPDATE : is this really necessary?
Individual::Individual(Individual &&other) noexcept
	: ID(std::move(other.ID)),
	  status(std::move(other.status)),
	  data(std::move(other.data)),
	  augm_map(std::move(other.augm_map)),
	  rate_map(std::move(other.rate_map)) {}

// deconstructor
Individual::~Individual() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Individual : implementations  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

//===== SETTER =====//

void Individual::assign_all_group_rates(const rateIDmapvec<const double *> &ptrs)
{
	for (const auto &pair : ptrs)
	{
		rateID r = pair.first;
		std::vector<const double *> vec_par_ptr = pair.second;

		if (vec_par_ptr.empty())
			std::invalid_argument(
				"Missing rate \"" + r + "\" for individual " + ID + ".");
		else
			this->rate_map[r] = GroupRate(r, vec_par_ptr);
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Individual : data  /////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if the individal has the data defined */
bool Individual::assert_data(const dataID &d) const
{
	if (this->data.find(d) == this->data.end())
		throw std::invalid_argument(
			"Unknown data \"" + d + "\" for individual " + ID + ".");
	else
		return true;

	return false;
}

/* Check if the given augmented variable is defined, otherwise raise error. */
bool Individual::assert_augm(const augmID &a) const
{
	if (this->augm_map.find(a) == this->augm_map.end())
		throw std::invalid_argument(
			"Unknown augmented variable \"" + a + "\" for individual " + ID + ".");
	else
		return true;

	return false;
}

/* Check if the given rate is defined, otherwise raise error. */
bool Individual::assert_rate(const rateID &r) const
{
	if (this->rate_map.find(r) == this->rate_map.end())
		throw std::invalid_argument(
			"Unknown group rate \"" + r + "\" for individual " + ID + ".");
	else
		return true;

	return false;
}
