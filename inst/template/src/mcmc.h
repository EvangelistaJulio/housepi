/*
 * mcmc.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef MCMC__H
#define MCMC__H

#include <string>
#include <vector>
#include <set>
#include <map>
#include <cmath>
#include <random>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>

#ifdef USE_RCPP
#include <Rcpp.h>
#endif

#include "individual.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  TYPEDEF  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

typedef std::map<string, string> ssargs;                

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
	//===== Setup =====//
	ulong setup_num_chains = 1;	  // n. chains
	ulong setup_num_iter = 10000; // n. iterations
	ulong setup_num_adapt = 0;	  // n. (finite) adaptive steps
	ulong setup_thinning = 1;	  // n. thinning iterations
	ulong setup_burn_in = 0;	  // n. burnin iterations

	bool setup_rnd_init_params = true; // random initialization of parameters
	unsigned long long setup_seed = 1; // random generator seed

	long double setup_inv_T = 1.0; // inverse temperature

	//===== Household Data =====//
	int num_indiv = 0;
	int num_househ = 0;
	vector<Household> househ = {};

	//----- Reports -----//
	std::map<string, int> report_data;
	std::map<string, std::map<string, int>> report_group_classes = {};

	//===== Augmented Variables =====//
	std::set<string> augm_full_set = {};
	std::set<string> augm_to_be_inferred = {};
	std::set<string> augm_with_model_set = {};
	std::set<string> augm_rand_init_set = {};

public:
	//===== CONSTRUCTOR =====//
	MCMC()
	{
		//----- default prng seed -----//
		auto time_now = std::chrono::steady_clock::now().time_since_epoch().count();
		setup_seed = static_cast<unsigned long long>(time_now);
	};
	virtual ~MCMC() {};

	//===== BUILD =====//
	void build_parameters(const string &pathtofile, bool verbose = true);
	void build_data(const string &pathtofile, bool verbose = true);
	void build_setup(const string &pathtofile, bool verbose = true);
	void build_augmented(const string &pathtofile, bool verbose = true);
	virtual void build_extra(const string &pathtofile, bool verbose = true);
	void build(const string &setup_file = UNDEFINED,
			   const string &params_file = UNDEFINED,
			   const string &data_file = UNDEFINED, 
			   const string &augm_file = UNDEFINED,
			   const string &extra_file = UNDEFINED, 
			   bool verbose = true);

	//===== SHOW =====//
	void show_setup();
	void show_summary();
	void show_group_classes();
	void show_parameters();
	void show_diagnostics();
	void show_results();

	//===== DATA =====//
	double get_data(int h, int i, const string d) const { return househ[h].indiv[i].get_data(d); }
	double get_data(int h, const string d) const { return househ[h].get_data(d); }

	//===== GROUP CLASSES =====//
	void set_rate_per_group(Individual &this_indiv);
	void set_rate_per_contact(int h);
	virtual void assign_rate_per_group_to_individual(Individual &this_indiv);
	virtual void assign_rate_per_contact_to_household(int h, int from, int to);
	double get_group_rate(int h, int i, const string r) const { return househ[h].indiv[i].get_group_rate(r); }
	double get_contact_rate(int h, int from, int to) const { return househ[h].get_contact_rate(from, to); }

	//===== AUGMENTED VARIABLES =====//
	void set_augm(int h, int i, const string a, double new_value) { househ[h].indiv[i].set_augm(a, new_value); };
	double get_augm(int h, int i, const string a) const { return househ[h].indiv[i].get_augm(a); }
	void assign_augmented_variables();
	void init_augmented_variables(bool is_init = true);
	virtual void update_first_case(int h);
	virtual void augm_store_dependecies(int h, int i);
	virtual void augm_update(int h, int i, string augmID, double new_val);
	virtual double augm_propose(int h, int i, string augmID);
	virtual double get_infec_time(int h, int i) const;

	//===== EVALUATE =====//
	void set_global_loglike(long double x) { global_loglike = x; };
	long double get_global_loglike() const { return global_loglike; };
	virtual double compute_house_loglike(int h);
	virtual long double compute_global_loglike();
	long double evaluate();

	//===== INFER =====//
	void disable_assert();
	virtual void init_model_for_inference() {};
	virtual void init_single_chain() {};
	void run_monte_carlo(bool progbar, ulong sweep = 0, bool is_init_augm = true);
	void infer(const string output_filename, bool progbar = true);

	//===== MARGINAL LOGLIKELIHOOD (POWER POSTERIORS) =====//
	double evaluate_marginal(int num_rungs, bool progbar = true, bool power_fractions = false);

	//===== SIMULATE =====//
	virtual void init_model_for_simulation() {};
	virtual void simulate_single_house(int h);
	void simulate();

	//----- Save augmented variables -----//
	void print_augm_variables(const string output_filename);

#ifdef USE_RCPP
	//===== RCPP EXPORTS =====//
	Rcpp::List export_setup();
	Rcpp::DataFrame export_params();
	Rcpp::DataFrame export_data();
#endif

private:
	//===== Global Loglikelihood =====//
	long double global_loglike = 0.0;

	//===== Internal =====//
	ulong _n_contacts_ = 0; // number of contacts

	ulong _chain_ = 0;			// current chain
	ulong _sample_size_ = 0;	// sample size
	ulong _mcmc_sweeps_ = 0;	// sweeps per chain
	ulong _tot_num_sweeps_ = 0; // total number of sweeps

	ulong _rung_ = 0;	   // current thermodynamic rung
	ulong _num_rungs_ = 0; // number of rungs

	bool _is_power_post_ = false;		// power posterior mode
	bool _is_data_augm_ = false;		// data augmentation mode
	bool _is_infer_complete_ = false;	// inference mode
	bool _is_print_inferstate_ = false; // print inference state mode

	vector<string> _data_headers_ = {};

	//===== Diagnostics =====//
	ulong dgn_ss = 0;			   // sample size
	long double dgn_avg = 0.;	   // sum of the loglikelihood
	long double dgn_sqr = 0.;	   // sum of the squared loglikelihood
	long double dgn_lppd = 0.;	   // log pointwise predictive density
	long double dgn_eff_pars = 0.; // effective n. of parameters

	long double dgn_var() { return (dgn_sqr - std::pow(dgn_avg, 2)) * dgn_ss / (dgn_ss - 1); }

	//===== BUILD =====//

	//----- build : auxiliary attributes -----//
	std::ofstream inferstate_ofstream;
	vector<ssargs> load_augm_definition = {};

	//----- build : auxiliary functions -----//
	bool stobool(const string &tkn);

	vector<string> parse_list_args(const string tkn, const char = screen::C_LIST_DELIM);
	dictargs merge_list_args(const string tkn_keys, const string tkn_vals);

	dictargs row_to_sdargs(ssargs &this_row);
	string string_cleaner(string tkn);
	vector<string> read_column_headers(std::ifstream &input,
												 const string pathtofile = UNDEFINED,
												 vector<string> required_headers = {});
	bool read_row_by_headers(std::ifstream &input, vector<string> headers, ssargs &this_row);
	void update_report_data();

	//----- build : type-specified variables -----//
	double type_to_double(const string tkn_val, const string tkn_typ, int h, int i);
	const double *type_to_doubleptr(const string tkn_val, const string tkn_typ, int h, int i);
	dictargs typedict_to_sdargs(const string tkn_key, const string tkn_val, const string tkn_typ, int h, int i);
	dictargs_ptr typedict_to_sdargs_ext(const string tkn_key, const string tkn_val, const string tkn_typ, int h, int i);

	//===== PRINTS =====//
	string get_progtxt(const string txt);

	//----- infer : save inference state -----//
	void inferstate_init(const string output_filename);
	void inferstate_print();
	void inferstate_close();

	//===== INFER : Monte Carlo Markov Chain =====//
	void init_move_counts();
	void init_diagnostics();
	void update_diagnostics();
	void close_diagnostics();

	void move_single_param(const string parID);
	void move_and_adapt_all_parameters(ulong step);
	void move_all_parameters();
	void move_all_augm_variables();
	long double move_single_house_augm(int h);

	double get_augm_accept_rate(const string augmID);
};

#endif // MCMC__H