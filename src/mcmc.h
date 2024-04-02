// mcmc.h //

#ifndef MCMC__H
#define MCMC__H

#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>

#include "individual.h"

// default acceptation rate for MCMC
#define DEFAULT_ACCEPT_RATE 0.24

// verbose state
extern bool Verbose;

//=====Forward declaration of the Infectivity Profile Function=====//
class Infectivity_Profile_Function;

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : Class ////////////////////////// //
////////////////////////////////////////////////////////////////////////

class MCMC : public Parameter_Dictionary
{
protected:
	bool setup_state = false;
	int setup_seed = 20210415;
	int setup_num_iter = 20000;
	int setup_num_adapt_iter = 5000;
	int setup_thinning = 20;
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> uniform;

	Individual *indiv;
	Household *househ;
	int num_indiv = 0;
	int num_househ = 0;

	const double **contact_patterns;

	Infectivity_Profile_Function *infec_profile_func;
	std::set<std::string> augm_set = {};
	std::set<std::string> augm_with_model_set = {};
	std::set<std::string> augm_to_infer_set = {};
	std::vector<std::string> data_headers = {};

public:
	//=====Constructors=====//
	MCMC() : uniform(0.0, 1.0) {}
	virtual ~MCMC(){};

	//=====deal with flatten matrices=====//
	int from_to(int i, int j) { return i * num_indiv + j; }
	double get_contact_from_to(int i, int k) { return (*contact_patterns[from_to(i, k)]); }

	//=====Build data and parameters=====//
	// load standard files
	void build_setup(std::string data_pathtofile);
	void build_parameters(std::string data_pathtofile);
	void build_augmented(std::string data_pathtofile);
	void build_data(std::string data_pathtofile);

	// load extra non standard files
	virtual void build_extra(std::string data_pathtofile);

	//=====Initialization=====//
	// attach individuals and assign group rate
	void add_individual(Individual new_indiv);
	int get_status(int i) const { return indiv[i].get_status(); }

	virtual double get_infec_time(int i) const;

	virtual void assign_rate_per_group_to_individual(int i);
	virtual void assign_rate_per_contact(int i, int k);
	virtual void assign_community_structure(std::string HH_index_column);

	virtual void assign_end_followup(int h);
	virtual void update_first_case(int h);

	void assign_augmented_variables();
	virtual void assign_infec_profile_function(){};
	virtual void initialize_model(){};

	//=====Log-Likelihood computations=====//
	virtual void store_augm_var_dependecies(int i);
	virtual double compute_single_house_loglike(int h);
	virtual double compute_global_loglike();
	void set_global_loglike(double x) { global_loglike = x; };
	double get_global_loglike() const { return global_loglike; };

	//=====Monte Carlo Markov Chain moves=====//
	void reset_moves();
	void move_single_param(std::string parID);
	void move_parameters();
	void move_and_adapt_parameters(int step, double target = DEFAULT_ACCEPT_RATE);
	void move_augm_variables();

	//=====infer/simulate=====//
	void inference(std::string output_filename, bool is_rand_params = false, double inv_temperature = 1.0);
	virtual void simulate_infection(){};
	void simulation(std::string output_filename);

	//=====Save parameters=====//
	void init_param_file(std::string output_filename);
	void print_parameters(int step);
	void close_param_file() { param_ofstream.close(); }

	//=====Save augmented variables=====//
	void print_augm_variables(std::string output_filename);

private:
	double global_loglike = 0.0;
	double inverse_temperature = 1.0;

	// outoput parameter file
	std::ofstream param_ofstream;
	std::map<std::string, std::map<std::string, int>> model_report = {};

	// accpetance rates
	std::unordered_map<std::string, int> m_num_accept_data = {};
	std::unordered_map<std::string, int> m_num_propos_data = {};

	std::vector<ssargs> load_augm_definition = {};

	//=====internal auxiliary functions=====//
	// file reading and formatting
	double type_to_double(std::string token, int i, char option_delim = OPTION_DELIM);
	const double *type_to_doubleptr(std::string token, char option_delim);
	sdargs typedict_to_sdargs(std::string token, int i, char dict_delim = DICT_DELIM, char option_delim = OPTION_DELIM);
	sdargs_ext typedict_to_sdargs_ext(std::string token, char dict_delim = DICT_DELIM, char option_delim = OPTION_DELIM);
};

#endif // MCMC__H
