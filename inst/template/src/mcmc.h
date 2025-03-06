/*
 * mcmc.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__MCMC__H
#define HOUSEPI__MCMC__H

#include <memory>
#include <string>
#include <vector>

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  MCMC : class  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Household models base class 'Monte Carlo Markov Chain'. */
class MCMC
{
public:
	//=======================//
	//===== CONSTRUCTOR =====//
	//=======================//

	MCMC();
	virtual ~MCMC();

	//================//
	//===== DATA =====//
	//================//

	Household &get_household(size_t H) const;
	size_t n_houses() const;
	size_t report_data(const std::string &key) const;

	//----- Group & Contact Rates -----//
	void build_all_rates();
	void emplace_group_rate(rateID r, parID p);
	void emplace_contact_rate(parID p);
	virtual void assign_group_rate_to_individual(Household &house, indiv i);
	virtual void assign_contact_rate_to_household(Household &house, indiv from, indiv to);
	void set_rate_per_group(Household &house, indiv i);
	void set_rate_per_contact(Household &house);

	//=================================//
	//===== RAND NUMBER GENERATOR =====//
	//=================================//

	double uniform();

	//======================//
	//===== PARAMETERS =====//
	//======================//

	double param(const parID &p) const;
	const double *param_pointer(const parID &p) const;

	//===============================//
	//===== AUGMENTED VARIABLES =====//
	//===============================//
	
	bool is_augm_any() const;
	bool is_augm_sample() const;
	bool is_augm_model() const;
	bool is_augm_generic() const;
	
	void init_augmented_variables();
	void assign_all_augm(Household &house, indiv i);
	void clear_all_augm(Household &house, indiv i);
	double compute_augm_model_loglike(const Household &house) const;

	//============================//
	//===== EXTRA ATTRIBUTES =====//
	//============================//

	void init_extra_dependecies();
	virtual void update_extra_dependecies(Household &house, indiv i);

	void allocate_extra_house_integer_vector(const extraID &e);
	void allocate_extra_house_double_vector(const extraID &e);
	void allocate_extra_house_double_matrix(const extraID &e);

	//===============================//
	//===== INFECTIVITY PROFILE =====//
	//===============================//

	virtual void init_infec_profile();
	void assign_default_infec_profile(const std::string &default_name);
	double infec_profile_inst(double gener_time, const dbl_vec &opts = {}) const;
	double infec_profile_cumul(double gener_time, const dbl_vec &opts = {}) const;

	//====================================//
	//===== EVALUATE (loglikelihood) =====//
	//====================================//

	virtual double compute_house_loglike(Household &house) const;
	void update_simple_global_loglike();
	void update_extra_global_loglike();
	double evaluate();

	//=================//
	//===== INFER =====//
	//=================//

	virtual void init_model_for_inference();
	void infer();

	//----- Marginal Loglikelihood (Power Posteriors) -----//
	double estimate_evidence(const std::string &scheme,
							 int n_rungs);

	//====================//
	//===== SIMULATE =====//
	//====================//

	virtual void init_model_for_simulation();
	virtual void simulate_single_house(Household &house);
	void simulate();

#ifdef _HOUSEPI_RCPP
	//========================//
	//===== Rcpp IMPORTS =====//
	//========================//

	void import_setup(const Rcpp::List &setup);
	void import_infectivity_profile(const Rcpp::List &setup);
	void import_parameter(const Rcpp::DataFrame &params);
	void import_dataset(const Rcpp::DataFrame &data);
	void import_augmented(const Rcpp::DataFrame &augm);

	//========================//
	//===== Rcpp EXPORTS =====//
	//========================//

	Rcpp::List export_setup() const;
	Rcpp::DataFrame export_parameter() const;
	Rcpp::DataFrame export_dataset() const;
	Rcpp::NumericMatrix export_traces() const;
	Rcpp::NumericMatrix export_household_loglike_matrix() const;
	Rcpp::DataFrame export_infectivity_profile(const Rcpp::List &opts,
											   const Rcpp::List &xlim,
											   int n_bins);

#else
	//================//
	//===== READ =====//
	//================//

	void read_setup(const std::string &pathtofile = NULL_STRING);
	void read_params(const std::string &pathtofile = NULL_STRING);
	void read_data(const std::string &pathtofile = NULL_STRING);
	void read_augm(const std::string &pathtofile = NULL_STRING);

	//================//
	//===== SHOW =====//
	//================//

	void show_summary() const;
	void show_full_build() const;
	void show_results() const;

	//----- C++ Print -----//
	void print_infer_state(const std::string &output_filename = NULL_STRING) const;
	void print_augm_variables(const std::string &output_filename = NULL_STRING) const;
#endif

private:
	//----- Pointer to Implementation -----//
	class Impl;
	std::unique_ptr<Impl> pImpl;

	//----- Error ----//
	std::string _error_message_ = "";

	//----- Group & Contact Rateas -----//
	rateIDmapvec<parID> _individual_group_rates_ = {};
	parID _pair_contact_rate_ = {};

	//----- Monte Carlo -----//
	void run_monte_carlo(progress::Bar &progbar, ulong sweep = 0);
	void move_all();
	void move_all_and_adapt_parameters(ulong step);
	bool move_independent_param(const parID &p);
	bool move_extra_depen_param(const parID &p);

	//----- Metropolis : augmented moves -----//
	double move_augm_single_house(Household &house);
	void metropolis_augm_single_move(Household &house, indiv i, const augmID &a, double &delta);

	void move_all_augm();

	void (MCMC::*_func_move_augm_indiv_)(Household &, indiv, double &) = &MCMC::_move_augm_indiv_generic_;
	void _move_augm_indiv_generic_(Household &house, indiv i, double &delta);
	void _move_augm_indiv_by_status_(Household &house, indiv i, double &delta);
};

#endif // HOUSEPI__MCMC__H
