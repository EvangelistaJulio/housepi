/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
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

	/* Virtual method for user initialization of aumented setup. */
	virtual void init_user_advanced_setup() {}

	//================//
	//===== DATA =====//
	//================//

	ulong get_report_indiv(const std::string &key) const;
	ulong get_report_house(const std::string &key) const;
	std::pair<ushort, ulong> get_report_hh_size(const std::string &key) const;

	//----- Group & Contact Rates -----//
	// UPDATE: group and contact rates should just be imported
	// i.e. not re-compilation
	void build_all_rates();

	/* Auxiliary wrapper to emplace_back */
	void emplace_group_rate(rateID r, parID p) { _individual_group_rates[r].emplace_back(p); }
	/* Auxiliary wrapper to emplace_back */
	void emplace_contact_rate(parID p) { _pair_contact_rate = p; }

	/* Assign group rate to an individual  */
	virtual void assign_group_rate_to_individual(Individual &) {}
	/* Assign the contact patterns between a couple of individuals. */
	virtual void assign_contact_rate_to_household(Household &, indiv, indiv) {}

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

	void assign_all_augm(Individual &member);
	void clear_all_augm(Individual &member);

	//============================//
	//===== CACHE ATTRIBUTES =====//
	//============================//

	/* Update cahced augmeted dependecies. */
	virtual void update_cache_dependecies(Household &, indiv) {}

	void assign_cache_house_integer_vector(const cacheID &e);
	void assign_cache_house_double_vector(const cacheID &e);
	void assign_cache_house_double_matrix(const cacheID &e);

	//===============================//
	//===== INFECTIVITY PROFILE =====//
	//===============================//

	void assign_template_infec_profile(const dict_str_str x);
	double infec_profile_inst(Individual &member, double t) const;
	double infec_profile_cumul(Individual &member, double t) const;

	//====================================//
	//===== EVALUATE (loglikelihood) =====//
	//====================================//

	/* Compute the loglikelihood for a single household. */
	virtual double compute_household_obs_loglike(Household &) const { throw std::runtime_error(
		"Called virtual loglikelihood method of undefined model."); }
	double evaluate();

	//=================//
	//===== INFER =====//
	//=================//

	/* Inference specific initialization (e. g. cache variables, lookup tables, etc). */
	virtual void init_model_for_inference() {}
	void infer();

	//====================//
	//===== SIMULATE =====//
	//====================//

	/* Simulation specific initialization (e. g. cache variables, lookup tables, etc). */
	virtual void init_model_for_simulation() {}

	/* Simulation for a single household. */
	virtual void simulate_single_household(Household &) { throw std::runtime_error(
		"Called virtual simulation method of undefined model."); }

#ifdef _HOUSEPI_RCPP

	//----- C++ Print -----//
	Rcpp::DataFrame simulate();
	Rcpp::List estimate_evidence(const std::string &scheme,
								 ulong n_rungs);

	//========================//
	//===== Rcpp IMPORTS =====//
	//========================//

	void import_setup(const Rcpp::List &input);
	void import_dataset(const Rcpp::DataFrame &input);
	void import_datalist(const Rcpp::DataFrame &input);
	void import_multiplicity(const Rcpp::DataFrame &input);
	void import_param(const Rcpp::DataFrame &input);
	void import_infec_profile(const Rcpp::List &input);
	void import_augm(const Rcpp::DataFrame &input);

	//========================//
	//===== Rcpp EXPORTS =====//
	//========================//

	Rcpp::List export_setup() const;
	Rcpp::DataFrame export_param() const;
	Rcpp::DataFrame export_dataset() const;
	Rcpp::DataFrame export_data_report() const;
	Rcpp::DataFrame export_house_size_report() const;
	Rcpp::NumericMatrix export_traces() const;
	Rcpp::NumericMatrix export_house_loglike_mtx() const;
	Rcpp::List showcase_infec_profile(int n_bins);

#else

	//----- C++ Print -----//
	void print_infer_state(const std::string &output_filename = HOUSEPI_NA_STRING) const;
	void print_simulated_dataset(const std::string &output_filename = HOUSEPI_NA_STRING);
	void print_evidence(const std::string &scheme, ulong n_rungs);

	//================//
	//===== READ =====//
	//================//

	void read_setup(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_dataset(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_datalist(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_multiplicity(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_param(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_infec_profile(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_augm(const std::string &pathtofile = HOUSEPI_NA_STRING);

	//================//
	//===== SHOW =====//
	//================//

	void show_setup() const;
	void show_all_params() const;
	void show_summary() const;
	void show_results() const;

#endif

private:
	//----- Pointer to Implementation -----//
	class Impl;
	std::unique_ptr<Impl> pImpl;

	//----- Group & Contact Rateas -----//
	rateIDmapvec<parID> _individual_group_rates = {};
	parID _pair_contact_rate = {};

	void _assign_group_rate(Individual &member);
	void _assign_contact_rate(Household &house);
};

#endif // HOUSEPI__MCMC__H
