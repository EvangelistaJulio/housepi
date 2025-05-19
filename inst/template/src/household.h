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

	std::string id = HOUSEPI_NA_STRING; // user id
	int status = MISSING_STATUS;		// status of the individual

	// NOTE: I am forced to store numerical variables
	dataIDmap<double> data_map = {};   // database associated
	augmIDmap<Variable> augm_map = {}; // augmented variables
	rateIDmap<GroupRate> rate_map;	   // group rates associated

	//=======================//
	//===== CONSTRUCTOR =====//
	//=======================//

	Individual(const std::string &id = HOUSEPI_NA_STRING,
			   int status = MISSING_STATUS,
			   const dataIDmap<double> &data = {});
	~Individual();

	//==================//
	//===== SETTER =====//
	//==================//

	void assign_all_group_rates(const rateIDmapvec<const double *> &ptrs);
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
	void assign_individual(const std::string &id,
						   int status,
						   const dataIDmap<double> &data);
	bool change_status(indiv i, int new_status);
	void assign_contact_rate(indiv from, indiv to, const double *x);

	//==================//
	//===== GETTER =====//
	//==================//

	const std::string &house_id() const;
	size_t index() const;
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

	//----- Individual : Setter -----//
	void assign_all_group_rates(indiv i, const rateIDmapvec<const double *> &x);

	//----- Individual : Getter -----//
	const std::string &id(indiv i) const;
	int status(indiv i) const;
	double data(indiv i, const dataID &d) const;
	double group_rate(indiv i, const rateID &r) const;
	double contact_rate(indiv from, indiv to) const;

	//=====================//
	//===== AUGMENTED =====//
	//=====================//

	//----- Augmented : Setter -----//
	void assign_augm(indiv i,
					 const augmID &a = HOUSEPI_NA_STRING,
					 const std::string &name = HOUSEPI_NA_STRING,
					 double value = HOUSEPI_NA_DOUBLE,
					 bool is_discrete = false,
					 const std::string &model_name = HOUSEPI_NA_STRING);

	void assign_augm_model_scalar(indiv i, const augmID &a,
								  const std::string &arg_name,
								  bool is_arg,
								  double arg_value);
	void assign_augm_model_pointer(indiv i, const augmID &a,
								   const std::string &arg_name,
								   bool is_arg,
								   const double *ptr_to_arg_value);

	void init_augm_model_reshaping(indiv i, const augmID &a);

	void assign_augm_proposer(indiv i,
							  const augmID &a = HOUSEPI_NA_STRING,
							  const std::string &rule_name = HOUSEPI_NA_STRING,
							  double rule_scale = HOUSEPI_NA_DOUBLE,
							  double rule_from = HOUSEPI_NA_DOUBLE);

	void set_augm_value(indiv i, const augmID &a, double x);
	void augm_model_rand(indiv i, const augmID &a, rnd_gen &prng);
	bool propose_valid_augm_move(indiv i, const augmID &a, rnd_gen &prng);
	void restore_augm_value(indiv i, const augmID &a);
	void clear_augm(indiv i, const augmID &a);

	//----- Augmented : Getter -----//
	double augm(indiv i, const augmID &a) const;
	const double *augm_pointer(indiv i, const augmID &a) const;
	bool is_augm_constant(indiv i, const augmID &a) const;
	bool is_augm_function(indiv i, const augmID &a) const;

	double augm_model_logpdf(indiv i, const augmID &a) const;
	double augm_model_logratio(indiv i, const augmID &a) const;
	double augm_propose_logratio(indiv i, const augmID &a) const;
	double augm_move_logratio(indiv i, const augmID &a) const;

	//============================//
	//===== EXTRA ATTRIBUTES =====//
	//============================//

	bool is_any_extra();
	void buffer_all_extra();
	void restore_all_extra();

	//----- Vector of integers -----//
	void allocate_extra_integer_vector(const extraID &e);
	int_vec &extra_integer_vector(const extraID &e);

	//----- Vector of doubles -----//
	void allocate_extra_double_vector(const extraID &e);
	dbl_vec &extra_double_vector(const extraID &e);

	//----- Matrix of doubles -----//
	void allocate_extra_double_matrix(const extraID &e);
	matrix<double> &extra_double_matrix(const extraID &e);

	//=========================//
	//===== LOGLIKELIHOOD =====//
	//=========================//

	//----- Loglikelihood : Setter -----//
	void set_augm_loglike(double x);
	void set_obs_loglike(double x);
	void increment_augm_loglike(double x);
	void restore_loglike();

	//----- Loglikelihood : Getter -----//
	double get_augm_loglike() const;
	double get_obs_loglike() const;
	double get_delta_obs_loglike() const;

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
	extraIDmap<int_vec> _extra_int_vec = {};
	extraIDmap<dbl_vec> _extra_dbl_vec = {};
	extraIDmap<matrix<double>> _extra_dbl_mtx = {};

	//----- Buffer EXTRA Attributes -----//
	extraIDmap<int_vec> _buffer_extra_int_vec = {};
	extraIDmap<dbl_vec> _buffer_extra_dbl_vec = {};
	extraIDmap<matrix<double>> _buffer_extra_dbl_mtx = {};
};

#endif // HOUSEPI__HOUSEHOLD__H
