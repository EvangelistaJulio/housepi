// mcmc.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

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

// to switch off this option, silent next line
#define PRINT_ADAPT_ITER // print also adaptive iterations

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////  MCMC : class  ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class MCMC : public Monte_Carlo_Core
{
protected:
	//=====Setup=====//
	int setup_seed = time(NULL);  // random seed
	int setup_num_iter = 20000;	  // number of iterations
	int setup_num_adapt_iter = 0; // number of (finite) adaptive steps
	int setup_thinning = 20;	  // number of thinning iterations

	//=====Database=====//
	std::vector<std::string> data_headers = {};

	//=====List of Households=====//
	std::vector<Household> househ = {};
	int num_househ = 0;

	//=====Augmented Variables=====//
	bool is_data_augmentation = false;
	std::set<std::string> augm_full_set = {};
	std::set<std::string> augm_to_be_inferred = {};
	std::set<std::string> augm_with_model_set = {};
	std::set<std::string> augm_rand_init_set = {};

	//=====Diagnostics=====//
	int dgn_ss = 0;						  // sample size
	long double dgn_sum_loglike = 0.;	  // sum of the loglikelihood
	long double dgn_sum_sqr_loglike = 0.; // sum of the squared loglikelihood

	//=====Reports (Database, Model)=====//
	std::map<std::string, int> report_data = {
		{"n.indiv", 0}, {"n.cases", 0}, {"n.unk.stat", 0}, {"no.HH.case", 0}, {"lst.1.case", 0}, {"lst.2.case", 0}, {"min.HH.size", 1000}, {"max.HH.size", 0}};
	std::map<std::string, std::map<std::string, int>> report_model = {};

public:
	//=====Build data and parameters=====//
	void build_setup(std::string data_pathtofile);
	void build_parameters(std::string data_pathtofile);
	void build_augmented(std::string data_pathtofile);
	void build_data(std::string data_pathtofile);	
	virtual void build_extra(std::string data_pathtofile); // load non standard file

	//=====General Initialization=====//
	virtual void assign_rate_per_group_to_individual(Individual &this_indiv);
	virtual void assign_rate_per_contact_to_household(int h, int from, int to);
	virtual void initialize_model() {};
	void set_rate_per_group(Individual &this_indiv);
	void set_rate_per_contact(int h);
	void display_report();

	//=====Augmented Variables=====//
	void assign_augmented_variables();
	void init_augmented_variables();
	virtual void update_first_case(int h);
	virtual void augm_var_store_dependecies(int h, int i);
	virtual double augm_var_propose(int h, int i, std::string augmID);
	virtual void augm_var_update(int h, int i, std::string augmID, double new_val);
	virtual double get_infec_time(int h, int i) const;

	//=====Log-Likelihood computations=====//
	void set_global_loglike(double x) { global_loglike = x; };
	double get_global_loglike() const { return global_loglike; };
	virtual double compute_single_house_loglike(int h);
	virtual double compute_global_loglike();

	//=====Monte Carlo Markov Chain moves=====//
	void reset_move_counts();
	void single_param_move(std::string parID);
	void move_all_parameters();
	double single_household_augm_move(int h);
	double get_augm_accept_rate(std::string augmID);
	void move_all_augm_variables();
	void move_and_adapt_all_parameters(int step);
	void update_model_diagnostics();
	void display_thin_results();

	//=====Infer=====//
	void inference(
		std::string output_filename,
		bool is_rand_params = false,
		double inv_temperature = 1.0);

	//=====Simulate=====//
	virtual void simulate_single_house(int h);
	virtual void format_augmented_data() {};
	void simulation();

	//=====Save parameters=====//
	void init_param_file(std::string output_filename);
	void print_parameters(int step);
	void close_param_file() { param_ofstream.close(); }

	//=====Save augmented variables=====//
	void print_augm_variables(std::string output_filename);

private:
	double global_loglike = 0.0;
	double inverse_temperature = 1.0;

	//=====BUILD=====//

	//-----build : auxiliary attributes-----//
	std::ofstream param_ofstream; 	// outoput parameter file
	std::vector<ssargs> load_augm_definition = {}; // temporary loader for augmented variable definition

	//----build : auxiliary functions-----//
	std::string string_cleaner(std::string token);
	bool stobool(std::string token);
	ssargs parse_dict_string(std::string token);
	sdargs parse_dict_args(std::string token);
	std::vector<std::string> read_column_headers(std::ifstream &input);
	bool read_row_by_headers(std::ifstream &input, std::vector<std::string> headers, ssargs &this_row);
	sdargs row_to_sdargs(ssargs &this_row);
	
	//----build : type-specified variables-----//
	double type_to_double(std::string token, int h, int i);
	const double *type_to_doubleptr(std::string token);
	sdargs typedict_to_sdargs(std::string token, int h, int i);
	sdargs_ext typedict_to_sdargs_ext(std::string token);
};

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

#endif // MCMC__H
