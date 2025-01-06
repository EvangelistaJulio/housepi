/*
 * individual.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__INDIVIDUAL__H
#define HOUSEPI__INDIVIDUAL__H

#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

#include "variable.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Group Rate : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

class GroupRate
{
public:
	std::string name = "";

	//===== CONSTRUCTOR ====//
	GroupRate(const std::string &name = UNDEFINED,
			  const std::vector<const double *> &pointers = {});
	~GroupRate();

	//===== GETTER ====//
	double get() const;

private:
	ushort n_rates = 0;
	const double ** pointers = {};
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Individual : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Struct for individual household study member. */
struct Individual
{
	//===== PUBLIC ATTRIBUTES =====//
	std::string ID = UNDEFINED;	 // user ID
	int status = 0;				 // status of the individual
	dataIDmap<double> data = {}; // database associated to individual

	rateIDmap<GroupRate> rate_map = {};
	augmIDmap<Base_Variable> augm_map = {};

	//===== CONSTRUCTOR =====//
	Individual(const std::string &ID,
			   int status,
			   const dataIDmap<double> &data = {});
	~Individual();

	void assign_all_group_rates(const rateIDmapvec<const double *> &group_rate_ptrs);

	//===== ASSERTION =====//
	// UPDATE : temporary since about to drop maps
	bool assert_data(const dataID &d) const;
	bool assert_augm(const augmID &a) const;
	bool assert_rate(const rateID &r) const;
};

#endif // HOUSEPI__INDIVIDUAL__H
