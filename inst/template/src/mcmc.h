/*
 * mcmc.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef MCMC__H
#define MCMC__H

#include <memory>
#include <map>
#include <cmath>
#include <random>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <iomanip>

#ifdef USE_RCPP
#include <Rcpp.h>
#endif

#include "household.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

// Forward declaration of implementation class
class Monte_Carlo_Build;

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  MCMC : class  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Household models base class 'Monte Carlo Markov Chain'. */
class MCMC
{
public:
	//===== CONSTRUCTOR =====//
	MCMC();
	~MCMC();

	//===== READ =====//
	void read_setup(const std::string &pathtofile = UNDEFINED);
	void read_params(const std::string &pathtofile = UNDEFINED);
	void read_data(const std::string &pathtofile = UNDEFINED);
	void read_augm(const std::string &pathtofile = UNDEFINED);

#ifdef USE_RCPP
	//===== Rcpp IMPORTS =====//
	void import_setup(Rcpp::List setup);
	void import_params(Rcpp::DataFrame params);
	void import_data(Rcpp::DataFrame data);
	void import_augm(Rcpp::DataFrame augm);

	//===== Rcpp EXPORTS =====//
	Rcpp::List export_setup() const;
	Rcpp::DataFrame export_params() const;
	Rcpp::DataFrame export_data() const;
	Rcpp::DataFrame export_traces() const;
	Rcpp::NumericMatrix export_loglike_matrix() const;
#else
	//===== Cpp SHOW =====//
	void show_full_build() const;
	void show_results() const;

	//===== Cpp PRINT =====//
	void print_infer_state(const std::string &output_filename = UNDEFINED) const;
	void print_augm_variables(const std::string &output_filename = UNDEFINED) const;
#endif

	//===== SHOW =====//
	void show_summary() const;
	void show_group_classes() const;

	//===== DATA =====//
	Household &house(size_t H);
	size_t get_num_house() const;
	int report_data(const std::string &key) const;

	// note : I didn't find a way to PIMPL these because of the virtual functions
	void build_data(const std::map<std::string, std::vector<std::string>> &datatable);
	void set_rate_per_group(Individual &this_indiv);
	void set_rate_per_contact(Household &this_house);

	virtual void assign_group_rate_to_individual(Individual &this_indiv);
	virtual void assign_contact_rate_to_household(Household &this_house, int from, int to);

	//===== RAND NUMBER GENERATOR =====//
	double uniform();

	//===== PARAMETERS =====//
	double get_param(const std::string &parID) const;
	const double *get_param_pointer(const std::string &parID) const;

	//===== AUGMENTED VARIABLES =====//
	void init_augmented_variables();
	virtual void update_augm_dependecies(Household &this_house, int i);
	void assign_augm(Individual &this_indiv);
	void clear_augm(Individual &this_indiv);

	//===== EVALUATE (loglikelihood) =====//
	virtual double compute_house_loglike(const Household &this_house) const;
	double evaluate();

	//===== INFER =====//
	virtual void init_model_for_inference() {};
	void infer(bool progbar = true);

	//===== MARGINAL LOGLIKELIHOOD (POWER POSTERIORS) =====//
	double evaluate_marginal(int num_rungs, bool adaptive = true, bool progbar = true);

	//===== SIMULATE =====//
	virtual void init_model_for_simulation() {};
	virtual void simulate_single_house(Household &this_house);
	void simulate(bool progbar = true);

private:
	//===== Pointer to Implementation =====//
	std::unique_ptr<Monte_Carlo_Build> pImpl;

	//===== Global Loglikelihood =====//
	double buffer = INF_DOUBLE;
	double global_loglike = 0.0;

	//===== Internal =====//
	ulong _chain_ = 0;			 // current chain
	ulong _tot_thin_sweeps_ = 0; // total number of sweeps
	ulong _rung_ = 0;			 // current thermodynamic rung
	ulong _num_rungs_ = 0;		 // number of rungs

	// errors
	bool _error_occurred_ = false;
	std::string _error_message_ = "";

	//===== PRINTS =====//
	std::string get_progtxt(const std::string &txt, bool is_power_post);

	//===== EVALUATE =====//
	double get_global_loglike() const { return global_loglike; };
	double get_delta_global_loglike() const { return global_loglike - buffer; }
	void restore_global_loglike() { global_loglike = buffer; }
	void increment_global_loglike(double delta) { global_loglike += delta; };
	void update_global_loglike();

	//===== INFER =====//

	//----- Monte Carlo -----//
	void run_monte_carlo(bool progbar, ulong sweep = 0, bool is_power_post = false);
	void move_single_param(const std::string &parID);
	void move_and_adapt_all_parameters(ulong step);
	void move_all_parameters();
	
	double move_single_house_augm(Household &this_house);
	void update_inver_temper(double inv_T);


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

	//----- Metropolis : augmented moves -----//
	void move_all_augm_variables() { (this->*_move_augm_method_)(); }
	void (MCMC::*_move_augm_method_)() = &MCMC::_move_no_augm_;
	void _move_no_augm_() {};
	void _move_all_augm_();
};

#endif // MCMC__H