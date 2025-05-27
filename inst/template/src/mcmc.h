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

	//================//
	//===== DATA =====//
	//================//

	Household &get_household(size_t H) const;
	size_t n_houses() const;
	ulong get_report_indiv(const std::string &key) const;
	ulong get_report_house(const std::string &key) const;
	std::pair<ushort, ulong> get_report_hh_size(const std::string &key) const;

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

	bool is_augm_model() const;
	void assign_all_augm(Household &house, indiv i);
	void clear_all_augm(Household &house, indiv i);

	//============================//
	//===== CACHE ATTRIBUTES =====//
	//============================//

	virtual void update_cache_dependecies(Household &house, indiv i);

	void allocate_cache_house_integer_vector(const cacheID &e);
	void allocate_cache_house_double_vector(const cacheID &e);
	void allocate_cache_house_double_matrix(const cacheID &e);

	//===============================//
	//===== INFECTIVITY PROFILE =====//
	//===============================//

	virtual void init_infec_profile();
	void assign_template_infec_profile(std::map<std::string, std::string> x);
	double infec_profile_inst(double gener_time) const;
	double infec_profile_cumul(double gener_time) const;

	//====================================//
	//===== EVALUATE (loglikelihood) =====//
	//====================================//

	virtual double compute_house_obs_loglike(Household &house) const;
	double update_house_loglike(size_t r, const bool is_augm_ll, const bool cache_depen);
	double evaluate();

	//=================//
	//===== INFER =====//
	//=================//

	virtual void init_model_for_inference();
	void infer();

	//----- Marginal Loglikelihood (Power Posteriors) -----//
	const std::vector<dbl_vec> power_posteriors(const std::string &scheme,
												ulong n_rungs);

	//====================//
	//===== SIMULATE =====//
	//====================//

	virtual void init_model_for_simulation();
	virtual void simulate_single_house(Household &house);

#ifdef _HOUSEPI_RCPP

	Rcpp::DataFrame simulate();
	Rcpp::List estimate_evidence(const std::string &scheme,
								 ulong n_rungs);

	//========================//
	//===== Rcpp IMPORTS =====//
	//========================//

	void import_setup(const Rcpp::List &input);
	void import_infec_profile(const Rcpp::List &input);
	void import_param(const Rcpp::DataFrame &input);
	void import_dataset(const Rcpp::DataFrame &input);
	void import_augm(const Rcpp::DataFrame &input);
	void import_multiplicity(const Rcpp::DataFrame &input);

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

	void print_evidence(const std::string &scheme, ulong n_rungs);

	//================//
	//===== READ =====//
	//================//

	void read_dataset(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_setup(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_param(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_infec_profile(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_augm(const std::string &pathtofile = HOUSEPI_NA_STRING);
	void read_multiplicity(const std::string &pathtofile = HOUSEPI_NA_STRING);

	//================//
	//===== SHOW =====//
	//================//

	void show_setup() const;
	void show_all_params() const;
	void show_summary() const;
	void show_results() const;

	//----- C++ Print -----//
	void print_infer_state(const std::string &output_filename = HOUSEPI_NA_STRING) const;
	void print_simulated_dataset(const std::string &output_filename = HOUSEPI_NA_STRING);
#endif

private:
	//----- Pointer to Implementation -----//
	class Impl;
	std::unique_ptr<Impl> pImpl;

	//----- Group & Contact Rateas -----//
	rateIDmapvec<parID> _individual_group_rates = {};
	parID _pair_contact_rate = {};

	//----- Simulation -----//
	void _aux_init_simulate();
};

#endif // HOUSEPI__MCMC__H
