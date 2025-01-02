/*
 * mcmc.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__MCMC__H
#define HOUSEPI__MCMC__H

#include <memory>
#include <map>
#include <cmath>
#include <random>
#include <unordered_set>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <functional>

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#endif

#include "variable.h"
#include "household.h"
#include "profiles.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  MCMC : class  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Forward declaration of build implementation. */
class Monte_Carlo_Build;

/* Forward declaration of infectivity profile. */
class Base_InfecProfile;

/* Household models base class 'Monte Carlo Markov Chain'. */
class MCMC
{
public:
	//===== CONSTRUCTOR =====//
	MCMC();
	virtual ~MCMC();

	//===== READ =====//
	void read_setup(const std::string &pathtofile = UNDEFINED);
	void read_params(const std::string &pathtofile = UNDEFINED);
	void read_data(const std::string &pathtofile = UNDEFINED);
	void read_augm(const std::string &pathtofile = UNDEFINED);

#ifdef _HOUSEPI_RCPP
	//===== Rcpp IMPORTS =====//
	void import_setup(const Rcpp::List &setup);
	void import_params(const Rcpp::DataFrame &params);
	void import_data(const Rcpp::DataFrame &data);
	void import_augm(const Rcpp::DataFrame &augm);

	//===== Rcpp EXPORTS =====//
	Rcpp::List export_setup() const;
	Rcpp::DataFrame export_params() const;
	Rcpp::DataFrame export_data() const;
	Rcpp::DataFrame export_traces() const;
	Rcpp::NumericMatrix export_household_loglike_traces() const;
#endif

	//===== SHOW =====//
	void show_full_build() const;
	void show_results() const;
	void show_summary() const;

	//===== C++ PRINT =====//
	void print_infer_state(const std::string &output_filename = UNDEFINED) const;
	void print_augm_variables(const std::string &output_filename = UNDEFINED) const;

	//===== DATA =====//
	size_t n_houses() const;
	size_t report_data(const std::string &key) const;

	// note : I didn't find a way to PIMPL these because of the virtual functions
	void build_data(const dict_str_vec &datatable);
	void set_rate_per_group(Household &house, indiv i);
	void set_rate_per_contact(Household &house);

	//==== GROUP & CONTACT RATES =====//
	void emplace_group_rate(rateID r, parID p);
	void emplace_contact_rate(parID p);
	virtual void assign_group_rate_to_individual(Household &house, indiv i);
	virtual void assign_contact_rate_to_household(Household &house, indiv from, indiv to);

	//===== RAND NUMBER GENERATOR =====//
	double uniform();

	//===== PARAMETERS =====//
	double param(const parID &p) const;
	const double *param_pointer(const parID &p) const;

	//===== AUGMENTED VARIABLES =====//
	void init_augmented_variables();
	virtual void update_augm_dependecies(Household &house, indiv i);
	void assign_all_augm(Household &house, indiv i);
	void clear_all_augm(Household &house, indiv i);

	//===== EXTRA ATTRIBUTES =====//
	void assign_extra_integer_to_all_houses(const std::vector<extraID> &names, int xinit = 0.0);
	void assign_extra_double_to_all_houses(const std::vector<extraID> &names, double xinit = 0.0);
	void assign_extra_vector_to_all_houses(const std::vector<extraID> &names, double xinit = 0.0);
	void assign_extra_matrix_to_all_houses(const std::vector<extraID> &names, double xinit = 0.0);

	//===== INFECTIVITY PROFILE =====//
	void assign_infectivity_profile(const std::string &name);
	double infec_profile_inst(const std::vector<double> &x);
	double infec_profile_cumul(const std::vector<double> &x);

	//===== EVALUATE (loglikelihood) =====//
	virtual double compute_house_loglike(Household &house) const;
	void update_global_loglike();
	double evaluate();

	//===== INFER =====//
	virtual void init_model_for_inference() {}
	void infer(bool is_progbar = true);

	//===== MARGINAL LOGLIKELIHOOD (POWER POSTERIORS) =====//
	double evaluate_marginal(int n_rungs, std::string scheme = "adaptive", bool is_progbar = true);

	//===== SIMULATE =====//
	virtual void init_model_for_simulation() {}
	virtual void simulate_single_house(Household &house);
	void simulate(bool is_progbar = true);

private:
	//===== POINTER TO IMPLEMENTATION =====//
	std::unique_ptr<Monte_Carlo_Build> pBuild;

	//===== GROUP & CONTACT RATES =====//
	rateIDmapvec<parID> _individual_group_rates_ = {};
	parID _pair_contact_rate_ = {};

	//===== INFECTIVITY PROFILE =====//
	std::unique_ptr<Base_InfecProfile> _infec_profile_;

	//===== PRINTS =====//
	const std::string fmt_progtxt(const std::string &txt, bool is_power_post = false);

	//===== INFER =====//
	ulong _chain_ = 1;			 // current chain
	ulong _tot_thin_sweeps_ = 0; // total number of sweeps
	ulong _rung_ = 0;			 // current thermodynamic rung
	ulong _n_rungs_ = 0;		 // number of rungs

	//----- Monte Carlo -----//
	void update_inver_temper(double inv_T);
	void run_monte_carlo(CustomProgressBar &progbar, ulong sweep = 0, bool is_power_post = false);
	
	//----- Metropolis : parameter moves -----//
	void move_all_parameters();
	void move_and_adapt_all_parameters(ulong step);
	void move_single_param(const parID &p);
	
	bool (MCMC::*_func_param_metropolis_)(const parID &) = &MCMC::_param_metropolis_1_;
	bool _param_metropolis_1_(const parID &p);
	bool _param_metropolis_0_(const parID &p);
	bool _param_metropolis_T_(const parID &p);

	//----- Metropolis : augmented moves -----//
	double move_augm_single_house(Household &house);
	void metropolis_augm_single_move(Household &house, indiv i, const augmID &a, double &delta);

	void (MCMC::*_func_move_all_augm_)() = &MCMC::_move_no_augm_;
	void _move_no_augm_() {};
	void _move_all_augm_();

	void (MCMC::*_func_move_augm_indiv_)(Household &, indiv, double &) = &MCMC::_move_augm_indiv_by_status_;
	void _move_augm_indiv_by_status_(Household &house, indiv i, double &delta);
	void _move_augm_indiv_generic_(Household &house, indiv i, double &delta);

	bool (MCMC::*_func_augm_metropolis_)(const Household &, indiv, const augmID &, double &) = &MCMC::_augm_metropolis_1_;
	bool _augm_metropolis_1_(const Household &house, indiv i, const augmID &a, double &delta);
	bool _augm_metropolis_0_(const Household &house, indiv i, const augmID &a, double &delta);
	bool _augm_metropolis_T_(const Household &house, indiv i, const augmID &a, double &delta);

	//===== ERROR =====//
	std::string _error_message_ = "";
};

#endif // HOUSEPI__MCMC__H
