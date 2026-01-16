/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__MCMC__H
#define HOUSEPI__MCMC__H

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

#include "tools/keys.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  *** MCMC : class ***  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Individual;
class Household;
class MCMCImpl;

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

	size_t get_max_hh_cases() const;
	void assert_no_missing_status() const;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : group & contact rates  //////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// UPDATE: also contact rates should just be imported
	void build_all_rates();
	/// Auxiliary wrapper to emplace_back */
	void emplace_contact_rate(std::string p) { _pair_contact_rate = p; }
	/// Assign the contact patterns between a couple of individuals.
	virtual void init_house_contact_rates(Household &, const Individual &, const Individual &) {}

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : parameters  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Return value of parameter.
	double par(const std::string &p) const;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : augmented variables  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void init_missing_augm(Individual &member);
	void clear_all_augm(Individual &member);
	void request_augm_discrete(const std::string &x);
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  MCMC : infectivity profile  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	void init_infec_profile_template(hsp::dict<std::string> x);
	double infec_profile_inst(const Individual &member, double t) const;
	double infec_profile_cumul(const Individual &member, double t) const;
	void request_infec_profile_discrete();

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : inference //////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Inference specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_inference() {}
	void infer();

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////  MCMC : evaluate log-likelihood  /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Update cahced augmeted dependecies.
	virtual void update_cache_dependecies(Household &, Individual &) {}
	/// Compute the log-likelihood for a single household.
	virtual double compute_house_obs_log_lik(Household &) { throw std::runtime_error(
		"Called undefined log-likelihood method.\nOverride it in your custom model implementation before proceeding."); }
	double evaluate_log_lik(const std::string &house_id = "");

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////  MCMC : simulation  ////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Simulation specific initialization (e. g. cache variables, lookup tables, etc).
	virtual void init_model_for_simulation() {}

	/// Simulation for a single household.
	virtual void simulate_single_house(Household &) { throw std::runtime_error(
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
	Rcpp::CharacterVector export_param_list_const() const;
	Rcpp::CharacterVector export_param_list_sample_standard() const;
	Rcpp::CharacterVector export_param_list_sample_cache() const;
	Rcpp::List export_infec_profile() const;
	Rcpp::CharacterVector export_augm_list_draws() const;
	Rcpp::CharacterVector export_augm_list_full() const;
	Rcpp::CharacterVector export_augm_list_export() const;
	Rcpp::CharacterVector export_augm_list_discrete() const;
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
	double log_evidence(const std::string &scheme, unsigned long n_rungs);

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

	/// Contact Rateas
	/// UPDATE: this should be dropped in favour of table definitions
	std::string _pair_contact_rate = {};

	void init_contact_rates_(Household &house);

#ifdef _HOUSEPI_RCPP
	Rcpp::DataFrame helper_export_data_wide_(const std::vector<std::string> &columns, bool is_status) const;
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
