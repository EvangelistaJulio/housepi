/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * household.h
 *
 */

#ifndef HOUSEPI__HOUSEHOLD__H
#define HOUSEPI__HOUSEHOLD__H

#include <memory>
#include <string>
#include <vector>

#include "tools.h"
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

/* Class for individual household study member. */
class Individual
{
public:
	//===== CONSTRUCTOR =====//
	Individual(size_t index = 0,
			   const std::string &id = HOUSEPI_NA_STRING,
			   int status = MISSING_STATUS,
			   const dataIDmap<double> &data = {});
	~Individual();

	//===== GETTER =====//
	const std::string &id() const;
	size_t index() const;
	int status() const;

	double data(const dataID &d) const;
	double group_rate(const rateID &r) const;
	double augm(const augmID &a) const;
	Variable &get_augm_variable(const augmID &a);

	double get_augm_model_loglike() const;

	//===== SETTER =====//
	void set_status(int s);
	void assign_all_group_rates(const rateIDmapvec<const double *> &ptrs);
	void clear_augm(const augmID &a);
	void assign_augm(const augmID &a,
					 const std::string &name,
					 double value,
					 bool is_discrete,
					 bool is_constant,
					 const std::string &func_name,
					 const std::string &model_name);
	void init_augm_model(const std::vector<augmID> &augm_model_list);

private:
	// individual group rate nested class
	struct GroupRate
	{
		std::string group_rate_name = "";
		ushort n_rates = 0;
		const double **pointers = {};

		GroupRate(const std::string &name = HOUSEPI_NA_STRING,
				  const std::vector<const double *> &pointers = {});
		~GroupRate();
		double get() const;
	};

	std::string _id = HOUSEPI_NA_STRING; // user id
	size_t _index = 0;					 // internal index
	int _status = MISSING_STATUS;		 // status of the individual

	// NOTE: I am forced to store numerical variables
	dataIDmap<double> _data_map = {};	 // database associated
	augmIDmap<Variable> _augm_map = {};	 // augmented variables
	rateIDmap<GroupRate> _rate_map = {}; // group rates associated

	std::vector<std::function<double()>> _augm_model_logpdfs = {};
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Household : class  ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Household class */
class Household
{
public:
	//=======================//
	//===== CONSTRUCTOR =====//
	//=======================//

	Household(const std::string &house_id = HOUSEPI_NA_STRING,
			  size_t index = SIZE_MAX);
	~Household();

	//==================//
	//===== SETTER =====//
	//==================//

	void set_multiplicity(ulong m, size_t class_idx);
	void assign_member(size_t index,
					   const std::string &id,
					   int status,
					   const dataIDmap<double> &data);
	bool change_status(indiv i, int new_status);
	void assign_contact_rate(indiv from, indiv to, const double *x);

	//==================//
	//===== GETTER =====//
	//==================//

	const std::string &house_id() const;
	size_t house_index() const;
	ushort size() const;

	ushort n_infec_status() const;
	ushort n_healthy_status() const;
	ushort n_missing_status() const;

	const std::vector<indiv> &infected_status_list() const;
	const std::vector<indiv> &healthy_status_list() const;
	const std::vector<indiv> &missing_status_list() const;

	ulong multiplicity() const;
	size_t multiplicity_class_index() const;
	bool is_multiplicity_head_class() const;

	//======================//
	//===== INDIVIDUAL =====//
	//======================//

	//----- Individual : Getter -----//
	const std::string &id(indiv i) const;
	size_t index(indiv i) const;
	int status(indiv i) const;
	double data(indiv i, const dataID &d) const;
	double group_rate(indiv i, const rateID &r) const;
	double contact_rate(indiv from, indiv to) const;
	Individual &get_member(indiv i);

	//----- Individual : Setter -----//
	void assign_all_group_rates(indiv i, const rateIDmapvec<const double *> &x);

	//============================//
	//===== EXTRA ATTRIBUTES =====//
	//============================//

	bool is_any_cache();
	void buffer_all_cache();
	void recover_all_cache();

	//----- Vector of integers -----//
	void allocate_cache_integer_vector(const cacheID &e);
	int_vec &cache_integer_vector(const cacheID &e);

	//----- Vector of doubles -----//
	void allocate_cache_double_vector(const cacheID &e);
	dbl_vec &cache_double_vector(const cacheID &e);

	//----- Matrix of doubles -----//
	void allocate_cache_double_matrix(const cacheID &e);
	matrix<double> &cache_double_matrix(const cacheID &e);

	//=====================//
	//===== AUGMENTED =====//
	//=====================//

	double augm(indiv i, const augmID &a) const;
	void set_augm_value(indiv i, const augmID &a, double x);
	
	//=========================//
	//===== LOGLIKELIHOOD =====//
	//=========================//

	//----- Loglikelihood : Setter -----//
	void set_augm_loglike(double x);
	void set_obs_loglike(double x);
	void recover_loglike_value();

	//----- Loglikelihood : Getter -----//
	double get_augm_loglike() const;
	double get_obs_loglike() const;

	//=======================//
	//===== DIAGNOSTICS =====//
	//=======================//

	void init_diagnostics(ulong sample_size);
	void update_diagnostics(ulong count);
	const dbl_vec &get_smpl_loglike() const;

private:
	std::string _house_id = HOUSEPI_NA_STRING; // user house id
	size_t _index = -1;						   // index of house
	ushort _house_size = 0;					   // n. of house members
	ulong _house_multiplicity = 1;			   // n. of equivalent houses
	size_t _multiplicity_class = -1;		   // index of multiplicity class

	std::vector<Individual> _members = {}; // members

	ushort _n_infec = 0;   // n. members with infected status
	ushort _n_healthy = 0; // n. memebers healthy status
	ushort _n_missing = 0; // n. members missing status

	// wheter members correspond to a given status or not
	std::vector<indiv> _infec_list = {};
	std::vector<indiv> _healthy_list = {};
	std::vector<indiv> _missing_list = {};

	// contact rate pointers
	matrix<const double *> _contact_rates;

	//===== LOGLIKELIHOOD =====//

	double _obs_loglike = 0;  // current observation loglikelihood
	double _augm_loglike = 0; // current augmented model loglikelihood

	double _obs_buffer = 0;	 // buffer observation loglikelihood
	double _augm_buffer = 0; // buffer augmented model loglikelihood

	dbl_vec _smpl_loglike = {}; // diagnostics

	//===== EXTRA ATTRIBUTES =====//

	//----- Extra Attributes -----//
	cacheIDmap<int_vec> _cache_int_vec = {};
	cacheIDmap<dbl_vec> _cache_dbl_vec = {};
	cacheIDmap<matrix<double>> _cache_dbl_mtx = {};

	//----- Buffer EXTRA Attributes -----//
	cacheIDmap<int_vec> _buffer_cache_int_vec = {};
	cacheIDmap<dbl_vec> _buffer_cache_dbl_vec = {};
	cacheIDmap<matrix<double>> _buffer_cache_dbl_mtx = {};
};

#endif // HOUSEPI__HOUSEHOLD__H
