/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/types.h"

/// UPDATE: can MCMC be redesigned to remove dependency on Rcpp?
#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  *** MCMC : class ***  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class AugmContainer;
class Individual;
class Household;
class MCMCImpl;

class MCMC
{
public:
	///////////////////////////
	//////  CONSTRUCTOR  //////
	///////////////////////////

	MCMC(const std::vector<std::string> &data_wide_keys,
		 const std::vector<std::string> &hh_size_keys,
		 const std::vector<std::string> &hh_size_vec_keys,
		 const std::vector<std::string> &hh_double_vec_keys,
		 const std::vector<std::string> &hh_double_mtx_keys,
		 const std::vector<std::string> &data_long_keys,
		 const std::vector<std::string> &param_keys,
		 const std::vector<std::string> &group_rate_keys,
		 const std::vector<std::string> &augm_keys);
	virtual ~MCMC();

	/// Virtual method for user initialization of advanced setup.
	virtual void init_advanced_setup() {}

	/// Uniform distributed in [0, 1) rand number generator.
	double uniform();

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////  MCMC : summaries  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	size_t get_max_hh_cases() const noexcept;
	/// Throw an exception if at least a member as missing status.
	void assert_no_missing_status() const;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : group & contact rates  //////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// UPDATE: also contact rates should just be imported
	void build_contact_rates();

	/// Assign the contact patterns between a couple of individuals.
	/// UPDATE: this will be dropped in favor of table definition.
	virtual std::string init_contact_rates(Household &, const Individual &, const Individual &) { return ""; }

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : parameters  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Return value of parameter.
	double par(size_t p) const;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : augmented variables  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void init_missing_augm(Individual &member);
	void request_discrete_augm(size_t a);
	void request_discrete_param(size_t p);

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : infectivity profile  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void init_infec_profile_template(hsp::dict<std::string> x);
	double infec_profile_inst(const AugmContainer &member, double t) const;
	double infec_profile_cumul(const AugmContainer &member, double t) const;
	void request_discrete_infec_profile();

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////  MCMC : evaluate log-likelihood  /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

    /// Return cached size_t value.
    size_t &cache_size(Household &house, size_t e);
    size_t cache_size(const Household &house, size_t e) const;

    /// Return cached size_t vector.
    std::vector<size_t> &cache_size_vector(Household &house, size_t e);
    const std::vector<size_t> &cache_size_vector(const Household &house, size_t e) const;

    /// Return cached double vector.
    std::vector<double> &cache_double_vector(Household &house, size_t e);
    const std::vector<double> &cache_double_vector(const Household &house, size_t e) const;

    /// Return cached double matrix.
    std::vector<std::vector<double>> &cache_double_matrix(Household &house, size_t e);
    const std::vector<std::vector<double>> &cache_double_matrix(const Household &house, size_t e) const;

	/// Update cahced augmeted dependecies.
	virtual void update_cache_dependecies(Household &, Individual &) {}
	/// Compute the log-likelihood for a single household.
	virtual double compute_house_obs_log_lik(Household &) = 0;

	/// Initialize model and evaluate log-likelihood.
	double evaluate_log_lik(std::string_view house_id = "");

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : inference //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Inference specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_evaluation() {}
	/// Inference specific initialization for each household.
	/// NOTE: it's run right after 'init_model_for_evaluation'.
	virtual void init_house_for_evaluation(Household &) {}
	void run_inference();


	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : simulation  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Simulation specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_simulation() {}

	/// Simulation for a single household.
	virtual void simulate_single_house(Household &) = 0;

#ifdef _HOUSEPI_RCPP
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////  MCMC : Rcpp operations  //////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	Rcpp::List log_evidence(const std::string &scheme, unsigned long n_rungs);

	Rcpp::DataFrame data_wide_for_simulation() const;
	Rcpp::List run_simulation(size_t n_rep) { return helper_simulate(n_rep, false); }
	Rcpp::List run_simulation_for_predict(size_t n) { return helper_simulate(n, true); }

	/// UPDATE: these do not need to be exposed
	Rcpp::List helper_simulate(size_t n, bool is_expec_measure);
	Rcpp::DataFrame helper_simulate_single_repetition();

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////  MCMC : Rcpp imports  ///////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void import_setup(SEXP input);
	void import_data_wide(SEXP input);
	void import_data_long(SEXP input);
	void import_multiplicity(SEXP input);
	void import_group_rate(SEXP input);
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
	Rcpp::List export_infec_profile() const;
	Rcpp::DataFrame export_augmented_data_wide() const;
	Rcpp::List export_house_snapshot(const std::string &house_id) const;
	Rcpp::NumericMatrix flush_total_log_lik_draws();
	Rcpp::NumericMatrix flush_param_draws();
	Rcpp::NumericMatrix flush_param_lprior_draws();
	Rcpp::NumericMatrix flush_augm_draws();
	Rcpp::NumericMatrix flush_house_log_lik_draws();

#else
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : C++ standalone print  ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void print_infer_state(const std::string &output_filename = "") const;
	void print_simulated_data_wide(const std::string &output_filename = "");
	double log_evidence(std::string_view scheme, unsigned long n_rungs);

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////  MCMC : C++ standalone read  ///////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void read_setup(const std::string &pathtofile = "");
	void read_data_wide(const std::string &pathtofile = "");
	void read_data_long(const std::string &pathtofile = "");
	void read_multiplicity(const std::string &pathtofile = "");
	void read_group_rate(const std::string &pathtofile = "");
	void read_param(const std::string &pathtofile = "");
	void read_infec_profile(const std::string &pathtofile = "");
	void read_augm(const std::string &pathtofile = "");
#endif

private:
	std::unique_ptr<MCMCImpl> pImpl_;
};
