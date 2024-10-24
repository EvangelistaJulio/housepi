/*
 * mcmc.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef MCMC__H
#define MCMC__H

#include <map>
#include <cmath>
#include <random>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>

#ifdef USE_RCPP
#include <Rcpp.h>
#endif

#include "read.h"
#include "individual.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  MCMC : class  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class MCMC : public Monte_Carlo_Core
{
protected:
	//===== Household Data =====//
	int num_indiv = 0;
	int num_house = 0;
	std::vector<Household> house = {};

	//----- Reports -----//
	std::map<std::string, std::string> report_pathtofile;
	std::map<std::string, int> report_data;
	std::map<std::string, std::map<std::string, int>> report_group_classes = {};

public:
	//===== CONSTRUCTOR =====//
	MCMC()
	{
		//----- default prng seed -----//
		auto time_now = std::chrono::steady_clock::now().time_since_epoch().count();
		setup_seed = static_cast<unsigned long long>(time_now);
	};

	//===== READ =====//
	void read_setup(const std::string &pathtofile = UNDEFINED);
	void read_params(const std::string &pathtofile = UNDEFINED);
	void read_data(const std::string &pathtofile = UNDEFINED);
	void read_augm(const std::string &pathtofile = UNDEFINED);

	//===== BUILD =====//
	void build_setup();
	void build_params();
	void build_data(const std::map<std::string, std::vector<std::string>> &datatable);
	void build_augm();

	//===== SHOW =====//
	void show_setup() const;
	void show_params() const;
	void show_summary() const;
	void show_group_classes() const;
	void show_full_build() const;

#ifdef USE_RCPP
	//----- Rcpp imports -----//
	void import_setup(Rcpp::List setup);
	void import_params(Rcpp::DataFrame params);
	void import_data(Rcpp::DataFrame data);
	void import_augm(Rcpp::DataFrame augm);

	//----- Rcpp exports -----//
	Rcpp::List export_setup() const;
	Rcpp::DataFrame export_params() const;
	Rcpp::DataFrame export_data() const;
	Rcpp::DataFrame export_traces() const;
	Rcpp::NumericMatrix export_loglike_matrix() const;

#else
	//----- cpp outputs -----//
	void show_results() const;
	void print_infer_state(const std::string &output_filename = UNDEFINED) const;
	void print_augm_variables(const std::string &output_filename = UNDEFINED) const;
#endif

	//===== GROUP CLASSES =====//
	void set_rate_per_group(Individual &new_indiv);
	void set_rate_per_contact(Household &this_house);
	virtual void assign_group_rate_to_individual(Individual &new_indiv);
	virtual void assign_contact_rate_to_household(Household &this_house, int from, int to);
	double get_group_rate(const Individual &this_indiv, const std::string &rateID) const { return this_indiv.get_group_rate(rateID); }
	double get_contact_rate(const Household &this_house, int from, int to) const { return this_house.get_contact_rate(from, to); }
	void update_report_data();

	//===== AUGMENTED VARIABLES =====//
	void init_augmented_variables();
	void assign_augm(Individual &this_indiv);
	void clear_augm(Individual &this_indiv);
	virtual void update_augm_dependecies(Household &this_house, int i);

	//===== EVALUATE =====//
	virtual double compute_house_loglike(const Household &this_house) const;
	double evaluate();

	//===== INFER =====//
	virtual void init_model_for_inference() {};

	void run_monte_carlo(bool progbar, ulong sweep = 0);
	void infer(bool progbar = true);
	void disable_assert();

	//===== MARGINAL LOGLIKELIHOOD (POWER POSTERIORS) =====//
	double evaluate_marginal(int num_rungs, bool adaptive = true, bool progbar = true);

	//===== SIMULATE =====//
	virtual void init_model_for_simulation() {};
	virtual void simulate_single_house(Household &this_house);
	void simulate(bool progbar = true);

private:
	//===== Global Loglikelihood =====//
	double buffer = INF_DOUBLE;
	double global_loglike = 0.0;

	//===== Diagnostics =====//
	ulong dgn_count = 0;
	std::vector<double> dgn_glob_loglike = {};

	//===== Setup =====//
	ulong setup_num_chains = 1;	  // n. chains
	ulong setup_num_iter = 10000; // n. iterations
	ulong setup_num_adapt = 0;	  // n. (finite) adaptive steps
	ulong setup_thinning = 1;	  // n. thinning iterations
	ulong setup_burn_in = 0;	  // n. burnin iterations

	bool setup_rnd_init_params = true; // random initialization of parameters
	unsigned long long setup_seed = 1; // random generator seed

	double setup_inv_T = 1.0; // inverse temperature

	//===== Augmented Variables =====//
	std::unordered_set<std::string> augm_full_set = {};
	std::unordered_set<std::string> augm_to_be_inferred = {};
	std::unordered_set<std::string> augm_with_model_set = {};

	//===== Internal =====//
	ulong _n_contacts_ = 0; // number of contacts

	ulong _chain_ = 0;			 // current chain
	ulong _sample_size_ = 0;	 // sample size
	ulong _thin_sweeps_ = 0;	 // thin sweeps per simulation
	ulong _tot_thin_sweeps_ = 0; // total number of sweeps

	ulong _rung_ = 0;	   // current thermodynamic rung
	ulong _num_rungs_ = 0; // number of rungs

	bool _is_undef_params_ = false; // undefined parameters valid only if infer
	bool _is_power_post_ = false;	// power posterior mode
	bool _is_move_augm_ = false;	// move augmented

	std::vector<std::string> _data_headers_ = {};

	// errors
	bool _error_occurred_ = false;
	std::string _error_message_ = "";

	//===== BUILD =====//

	//----- build : auxiliary attributes -----//
	std::map<std::string, std::string> _setup_deflist_ = {};
	std::map<std::string, std::vector<std::string>> _param_deftable_ = {};
	std::map<std::string, std::vector<std::string>> _augm_deftable_ = {};

	//----- build : type-specified variables -----//
	double type_to_double(const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv);
	const double *type_to_doubleptr(const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv);
	dictargs typedict_to_sdargs(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv);
	dictargs_ptr typedict_to_sdargs_ptr(const std::string tkn_key, const std::string tkn_val, const std::string tkn_typ, const Individual &this_indiv);

	//===== PRINTS =====//
	std::string get_progtxt(const std::string &txt);

	//===== EVALUATE =====//
	void update_house_loglike(Household &this_house) { this_house.set_loglike(compute_house_loglike(this_house)); }
	double compute_augm_loglike(const Household &this_house) const;
	double get_global_loglike() const { return global_loglike; };
	double get_delta_global_loglike() const { return global_loglike - buffer; }
	void restore_global_loglike() { global_loglike = buffer; }
	void increment_global_loglike(double delta) { global_loglike += delta; };
	void update_global_loglike();

	//===== INFER =====//

	//----- Monte Carlo -----//
	void move_single_param(const std::string &parID);
	void move_and_adapt_all_parameters(ulong step);
	void move_all_parameters();
	void move_all_augm_variables();
	double move_single_house_augm(Household &this_house);
	void set_inver_temper(double inv_T);

	//----- Metropolis : parameter moves -----//
	bool metropolis_param_move(const std::string &parID) { return (this->*_param_metropolis_method_)(parID); }
	bool (MCMC::*_param_metropolis_method_)(const std::string &) = &MCMC::_param_metropolis_1_;
	bool _param_metropolis_1_(const std::string &parID);
	bool _param_metropolis_0_(const std::string &parID);
	bool _param_metropolis_T_(const std::string &parID);

	//----- Metropolis : augmented moves -----//
	bool metropolis_augm_move(const Household &this_house, int i, const std::string &augmID, double &delta)
	{
		return (this->*_augm_metropolis_method_)(this_house, i, augmID, delta);
	}
	bool (MCMC::*_augm_metropolis_method_)(const Household &, int, const std::string &, double &) = &MCMC::_augm_metropolis_1_;
	bool _augm_metropolis_1_(const Household &this_house, int i, const std::string &augmID, double &delta);
	bool _augm_metropolis_0_(const Household &this_house, int i, const std::string &augmID, double &delta);
	bool _augm_metropolis_T_(const Household &this_house, int i, const std::string &augmID, double &delta);

	//----- Values for Diagnostics and Results -----//
	void init_diagnostics();
	void update_diagnostics();
	void init_move_counts();
};

#endif // MCMC__H