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

#include "variable.h"

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
class Individual
{
public:
	//===== CONSTRUCTOR =====//
	Individual(const std::string &ID = NULL_STRING,
			   int status = MISSING_STATUS,
			   const dataIDmap<double> &data = {});
	Individual(Individual &&other) noexcept; // Move constructor
	~Individual();

	//===== PUBLIC ATTRIBUTES =====//
	std::string ID = NULL_STRING; // user ID
	int status = MISSING_STATUS;  // status of the individual

	//===== DATA =====//
	// UPDATE: if I had the chance to use boost I would do:
	// NOTE: BH R library doesn't contain variant
	// > #include <boost/variant.hpp>
	// > ...
	// > using ValueType = boost::variant<double, std::string>;
	// This way I could have string data (nice group_rates names)
	// The the getter would become something like:
	// template<typename T>
	// T get(const dataID& d) const { return boost::get<T>(data.find(d)->second); }
	//
	dataIDmap<double> data = {}; // database associated

	//===== AUGMENTED VARIABLES =====//
	augmIDmap<Base_Variable> augm_map = {}; // augmented variables

	//===== GROUP RATES =====//
	struct GroupRate
	{
		std::string name = "";
		ushort n_rates = 0;
		const double **pointers = {};
		
		GroupRate(const std::string &name = NULL_STRING,
				  const std::vector<const double *> &pointers = {});
		~GroupRate();
		double get() const;
	};

	rateIDmap<GroupRate> rate_map; // group rates associated

	void assign_all_group_rates(const rateIDmapvec<const double *> &ptrs);

	//===== ASSERTION =====//
	// UPDATE : assertion is a drawback, how about dropping maps?
	bool assert_data(const dataID &d) const;
	bool assert_augm(const augmID &a) const;
	bool assert_rate(const rateID &r) const;
};

#endif // HOUSEPI__INDIVIDUAL__H
