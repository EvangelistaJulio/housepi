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

	size_t n_houses() const;
	size_t report_data(const std::string &key) const;

	// note : I didn't find a way to PIMPL build_data
	// because of the virtual functions set_rate_per_group and set_rate_per_contact
	void build_data(const dict_str_vec &datatable);

	//----- Group & Contact Rates -----//
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

	void init_augmented_variables();
	virtual void update_augm_dependecies(Household &house, indiv i);
	void assign_all_augm(Household &house, indiv i);
	void clear_all_augm(Household &house, indiv i);

	//============================//
	//===== EXTRA ATTRIBUTES =====//
	//============================//

	void assign_extra_integer_to_all_houses(const std::vector<extraID> &names, int xinit = 0.0);
	void assign_extra_double_to_all_houses(const std::vector<extraID> &names, double xinit = 0.0);
	void assign_extra_vector_to_all_houses(const std::vector<extraID> &names, double xinit = 0.0);
	void assign_extra_matrix_to_all_houses(const std::vector<extraID> &names, double xinit = 0.0);

	//===============================//
	//===== INFECTIVITY PROFILE =====//
	//===============================//

	void assign_infec_profile(const std::string &default_name);
	double infec_profile_inst(double gener_time, const dbl_vec &opts = {}) const;
	double infec_profile_cumul(double gener_time, const dbl_vec &opts = {}) const;

	//====================================//
	//===== EVALUATE (loglikelihood) =====//
	//====================================//

	virtual double compute_house_loglike(Household &house) const;
	void update_global_loglike();
	double evaluate();

	//=================//
	//===== INFER =====//
	//=================//

	virtual void init_model_for_inference();
	void infer(bool is_progbar = true);

	//----- Marginal Loglikelihood (Power Posteriors) -----//
	double evaluate_marginal(int n_rungs, std::string scheme = "adaptive", bool is_progbar = true);

	//====================//
	//===== SIMULATE =====//
	//====================//

	virtual void init_model_for_simulation();
	virtual void simulate_single_house(Household &house);
	void simulate(bool is_progbar = true);

	//================//
	//===== SHOW =====//
	//================//

	void show_full_build() const;
	void show_results() const;
	void show_summary() const;

	//----- C++ Print -----//
	void print_infer_state(const std::string &output_filename = NULL_STRING) const;
	void print_augm_variables(const std::string &output_filename = NULL_STRING) const;

	//================//
	//===== READ =====//
	//================//

	void read_setup(const std::string &pathtofile = NULL_STRING);
	void read_params(const std::string &pathtofile = NULL_STRING);
	void read_data(const std::string &pathtofile = NULL_STRING);
	void read_augm(const std::string &pathtofile = NULL_STRING);

#ifdef _HOUSEPI_RCPP
	//========================//
	//===== Rcpp IMPORTS =====//
	//========================//

	void import_setup(const Rcpp::List &setup);
	void import_params(const Rcpp::DataFrame &params);
	void import_data(const Rcpp::DataFrame &data);
	void import_augm(const Rcpp::DataFrame &augm);

	//========================//
	//===== Rcpp EXPORTS =====//
	//========================//

	Rcpp::List export_setup() const;
	Rcpp::DataFrame export_params() const;
	Rcpp::DataFrame export_data() const;
	Rcpp::DataFrame export_traces() const;
	Rcpp::NumericMatrix export_household_loglike_traces() const;
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
	void run_monte_carlo(CustomProgressBar &progbar, ulong sweep = 0);
	void move_all();
	void move_all_and_adapt_parameters(ulong step);
	void move_single_param(const parID &p);

	//----- Metropolis : augmented moves -----//
	double move_augm_single_house(Household &house);
	void metropolis_augm_single_move(Household &house, indiv i, const augmID &a, double &delta);

	void (MCMC::*_func_move_all_augm_)() = &MCMC::_move_no_augm_;
	void _move_no_augm_();
	void _move_all_augm_();

	void (MCMC::*_func_move_augm_indiv_)(Household &, indiv, double &) = &MCMC::_move_augm_indiv_generic_;
	void _move_augm_indiv_generic_(Household &house, indiv i, double &delta);
	void _move_augm_indiv_by_status_(Household &house, indiv i, double &delta);
};

#endif // HOUSEPI__MCMC__H
