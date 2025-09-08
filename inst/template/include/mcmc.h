/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc.h
 *
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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  MCMC : class  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Household models base class 'Monte Carlo Markov Chain'.
class MCMC
{
public:
	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : constructor  ///////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	MCMC();
	virtual ~MCMC();

	/// Virtual method for user initialization of advanced setup.
	virtual void init_advanced_setup() {}

	/// Uniform distributed in [0, 1) rand number generator.
	double uniform();

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////  MCMC : summaries  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	unsigned long get_report_indiv(const std::string &s) const;
	unsigned long get_report_house(const std::string &s) const;
	std::pair<unsigned short, unsigned long> get_report_hh_size(const std::string &s) const;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : group & contact rates  //////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// UPDATE: also contact rates should just be imported
	void build_all_rates();
	/// Auxiliary wrapper to emplace_back */
	void emplace_contact_rate(std::string p) { _pair_contact_rate = p; }
	/// Assign the contact patterns between a couple of individuals.
	virtual void init_household_contact_rates(Household &, const Individual &, const Individual &) {}

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : parameters  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double param(const std::string &p) const;
	observer_dbl_ptr param_pointer(const std::string &p) const;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : augmented variables  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void init_all_augm(Individual &member);
	void clear_all_augm(Individual &member);

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////  MCMC : cache attributes  //////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Update cahced augmeted dependecies.
	virtual void update_cache_dependecies(Household &, Individual &) {}

	void cache_house_integer_vector(const std::string &e);
	void cache_house_size_t_vector(const std::string &e);
	void cache_house_double_vector(const std::string &e);
	void cache_house_double_matrix(const std::string &e);

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : infectivity profile  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////


	void init_infec_profile_template(dict_str x);
	double infec_profile_inst(const Individual &member, double t) const;
	double infec_profile_cumul(const Individual &member, double t) const;

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : inference //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Inference specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_inference() {}
	void infer();
	void clear_memory();

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////  MCMC : evaluate loglikelihood  //////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Compute the loglikelihood for a single household.
	virtual double compute_household_obs_loglike(Household &) { throw std::runtime_error(
		"Called virtual loglikelihood method of undefined model."); }
	double evaluate();

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : simulation  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Simulation specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_simulation() {}

	/// Simulation for a single household.
	virtual void simulate_single_household(Household &) { throw std::runtime_error(
		"Called virtual simulation method of undefined model."); }

#ifdef _HOUSEPI_RCPP
	//----- C++ Print -----//
	Rcpp::DataFrame data_wide_for_simulation() const;
	Rcpp::List simulate(size_t n_rep);
	Rcpp::List evidence(const std::string &scheme, unsigned long n_rungs);

	//========================//
	//===== Rcpp IMPORTS =====//
	//========================//

	void import_setup(SEXP input);
	void import_data_wide(SEXP input);
	void import_data_long(SEXP input);
	void import_multiplicity(SEXP input);
	void import_group_rates(SEXP input);
	void import_param(SEXP input);
	void import_infec_profile(SEXP input);
	void import_augm(SEXP input);

	//========================//
	//===== Rcpp EXPORTS =====//
	//========================//

	Rcpp::DataFrame export_data_wide() const;
	Rcpp::DataFrame export_data_long() const;
	Rcpp::DataFrame summarise_cases() const;
	Rcpp::DataFrame summarise_house_sizes() const;
	Rcpp::List export_setup() const;
	Rcpp::List export_infec_profile() const;
	Rcpp::DataFrame export_param() const;
	Rcpp::DataFrame export_augmented_data_wide() const;
	Rcpp::NumericMatrix export_draws() const;
	Rcpp::NumericMatrix export_house_loglike_draws() const;
	Rcpp::List export_house_snapshot(const std::string &house_id) const;
#else
	//----- C++ Print -----//
	void print_infer_state(const std::string &output_filename = HSP_NA_STRING) const;
	void print_simulated_data_wide(const std::string &output_filename = HSP_NA_STRING);
	void print_evidence(const std::string &scheme, unsigned long n_rungs);

	//////////////////////////
	/////////  read  /////////
	//////////////////////////

	void read_setup(const std::string &pathtofile = HSP_NA_STRING);
	void read_data_wide(const std::string &pathtofile = HSP_NA_STRING);
	void read_data_long(const std::string &pathtofile = HSP_NA_STRING);
	void read_multiplicity(const std::string &pathtofile = HSP_NA_STRING);
	void read_group_rates(const std::string &pathtofile = HSP_NA_STRING);
	void read_param(const std::string &pathtofile = HSP_NA_STRING);
	void read_infec_profile(const std::string &pathtofile = HSP_NA_STRING);
	void read_augm(const std::string &pathtofile = HSP_NA_STRING);

	//////////////////////////
	/////////  show  /////////
	//////////////////////////

	void show_setup() const;
	void show_all_params() const;
	void show_summary() const;
	void show_results() const;

#endif

private:
	class Impl; 	/// pointer to Implementation
	std::unique_ptr<Impl> pImpl_;

	/// Contact Rateas
	/// UPDATE: this should be dropped in favour of table definitions
	std::string _pair_contact_rate = {};

	void init_contact_rates_(Household &house);
};

#endif // HOUSEPI__MCMC__H
