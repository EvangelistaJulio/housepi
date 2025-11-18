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

#include "tools/types.h"
#include "tools/na_utils.h"
#include "household.h"

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  MCMC : class  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Report for all individuals in the data.
struct ReportIndividual
{
	size_t total = 0;
	size_t missing = 0;
	size_t healthy = 0;
	size_t infected = 0;
};

/// Report for all households in the data.
struct ReportHousehold
{
	size_t total = 0;
	size_t no_cases = 0;
	size_t at_least_1_case = 0;
	size_t at_least_2_cases = 0;
};

/// Report for all household sizes.
struct ReportHouseSizes
{
	std::pair<unsigned short, size_t> max_cases = {0, 0};
	std::pair<unsigned short, size_t> max_miss = {0, 0};
	std::pair<unsigned short, size_t> min_size = {1000, 0};
	std::pair<unsigned short, size_t> max_size = {0, 0};
};

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

	const ReportIndividual &get_report_indiv() const;
	const ReportHousehold &get_report_house() const;
	const ReportHouseSizes &get_report_hh_size() const;
	void assert_no_missing_status() const;

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

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : augmented variables  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void init_missing_augm(Individual &member);
	void clear_all_augm(Individual &member);
	void request_augm_discrete(const std::string &x);

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

	void init_infec_profile_template(u_dict<std::string> x);
	double infec_profile_inst(const Individual &member, double t) const;
	double infec_profile_cumul(const Individual &member, double t) const;
	void request_infec_profile_discrete();

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : inference //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Inference specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_inference() {}
	void infer();
	void clear_memory();

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////  MCMC : evaluate log-likelihood  /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Compute the log-likelihood for a single household.
	virtual double compute_household_obs_log_lik(Household &) { throw std::runtime_error(
		"Called undefined log-likelihood method.\nOverride it in your custom model implementation before proceeding."); }
	double evaluate_log_lik(const std::string &house_id = hsp::na::STRING);

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : simulation  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Simulation specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_simulation() {}

	/// Simulation for a single household.
	virtual void simulate_single_household(Household &) { throw std::runtime_error(
		"Called undefined simulation method.\nOverride it in your custom model implementation before proceeding."); }

#ifdef _HOUSEPI_RCPP
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////  MCMC : Rcpp operations  //////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	Rcpp::DataFrame data_wide_for_simulation() const;
	Rcpp::List simulate(size_t n_rep) { return helper_simulate_(n_rep, false); }
	Rcpp::List simulate_for_prediction(size_t n) { return helper_simulate_(n, true); }
	Rcpp::List log_evidence(const std::string &scheme, unsigned long n_rungs);

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  MCMC : Rcpp imports  ///////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void import_setup(SEXP input);
	void import_data_wide(SEXP input);
	void import_data_long(SEXP input);
	void import_multiplicity(SEXP input);
	void import_group_rates(SEXP input);
	void import_param(SEXP input);
	void import_infec_profile(SEXP input);
	void import_augm(SEXP input);
	void import_param_for_prediction(SEXP input);

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  MCMC : Rcpp exports  ///////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	Rcpp::DataFrame export_data_wide() const;
	Rcpp::DataFrame export_data_long() const;
	Rcpp::DataFrame summarise_cases() const;
	Rcpp::DataFrame summarise_house_sizes() const;
	Rcpp::List export_setup() const;
	Rcpp::CharacterVector export_param_list_const() const;
	Rcpp::CharacterVector export_param_list_sample_standard() const;
	Rcpp::CharacterVector export_param_list_sample_cache() const;
	Rcpp::List export_infec_profile() const;
	Rcpp::CharacterVector export_augm_list_draws() const;
	Rcpp::CharacterVector export_augm_list_full() const;
	Rcpp::CharacterVector export_augm_list_export() const;
	Rcpp::CharacterVector export_augm_list_discrete() const;
	Rcpp::DataFrame export_augmented_data_wide() const;
	Rcpp::NumericMatrix export_draws() const;
	Rcpp::NumericMatrix export_house_log_lik_draws() const;
	Rcpp::List export_house_snapshot(const std::string &house_id) const;
#else
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : C++ standalone print  ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void print_infer_state(const std::string &output_filename = hsp::na::STRING) const;
	void print_simulated_data_wide(const std::string &output_filename = hsp::na::STRING);
	double log_evidence(const std::string &scheme, unsigned long n_rungs);

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////  MCMC : C++ standalone read  ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void read_setup(const std::string &pathtofile = hsp::na::STRING);
	void read_data_wide(const std::string &pathtofile = hsp::na::STRING);
	void read_data_long(const std::string &pathtofile = hsp::na::STRING);
	void read_multiplicity(const std::string &pathtofile = hsp::na::STRING);
	void read_group_rates(const std::string &pathtofile = hsp::na::STRING);
	void read_param(const std::string &pathtofile = hsp::na::STRING);
	void read_infec_profile(const std::string &pathtofile = hsp::na::STRING);
	void read_augm(const std::string &pathtofile = hsp::na::STRING);
#endif

private:
	class Impl; /// pointer to Implementation
	std::unique_ptr<Impl> pImpl_;

	/// Contact Rateas
	/// UPDATE: this should be dropped in favour of table definitions
	std::string _pair_contact_rate = {};

	void init_contact_rates_(Household &house);

#ifdef _HOUSEPI_RCPP
	Rcpp::DataFrame helper_export_data_wide_(const str_vec &columns, bool is_status) const;
	Rcpp::NumericVector helper_param_rcpp_snapshot_() const;
	Rcpp::DataFrame helper_data_wide_rcpp_snapshot_(const Household &house) const;
	Rcpp::DataFrame helper_data_long_rcpp_snapshot_(const Household &house) const;
	Rcpp::DataFrame helper_group_rates_rcpp_snapshot_(const Household &house) const;
	Rcpp::List helper_infec_profile_rcpp_snapshot_(const Household &house) const;
	Rcpp::DataFrame helper_augm_rcpp_snapshot_(const Household &house) const;

	Rcpp::List helper_simulate_(size_t n, bool is_expec_measure);
	Rcpp::DataFrame helper_simulate_single_repetition_();
#endif
};

#endif // HOUSEPI__MCMC__H
