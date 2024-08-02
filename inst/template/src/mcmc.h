// mcmc.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#ifndef MCMC__H
#define MCMC__H

#include <string>
#include <vector>
#include <set>
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
	ulong setup_num_chains = 1;	  // n. chains
	ulong setup_num_iter = 10000; // n. iterations
	ulong setup_num_adapt = 0;	  // n. (finite) adaptive steps
	ulong setup_thinning = 1;	  // n. thinning iterations
	ulong setup_burn_in = 0;	  // n. burnin iterations
	ulong setup_sample_size = 0;  // n. prints (just internal)

	unsigned long long setup_seed = 1;	// random generator seed
	double setup_inv_temper = 1.0;		// inverse temperature
	bool setup_rnd_init_params = false; // random initialization of parameters

	//=====Internal=====//
	ulong intern_chain = 0;
	bool intern_is_infer_complete = false;
	bool intern_is_diagnostics = false;

	//=====Database=====//
	std::vector<std::string> data_headers = {};

	//=====List of Households=====//
	std::vector<Household> househ = {};
	int num_househ = 0;

	//=====Augmented Variables=====//
	bool setup_data_augm = false; // data augmentation mode
	std::set<std::string> augm_full_set = {};
	std::set<std::string> augm_to_be_inferred = {};
	std::set<std::string> augm_with_model_set = {};
	std::set<std::string> augm_rand_init_set = {};

	//=====Diagnostics=====//
	long double dgn_loglike = 0.;	  // sum of the loglikelihood
	long double dgn_sqr_loglike = 0.; // sum of the squared loglikelihood

	//=====Reports (Database, Model)=====//
	std::map<std::string, int> report_data;
	std::map<std::string, std::map<std::string, int>> report_model = {};

public:
	MCMC()
	{
		// default prng seed
		auto time_now = std::chrono::steady_clock::now().time_since_epoch().count();
		setup_seed = static_cast<unsigned long long>(time_now);
	};
	virtual ~MCMC() {};

	//=====Build data and parameters=====//
	void build_parameters(std::string pathtofile, bool verbose = true);
	void build_data(std::string pathtofile, bool verbose = true);
	void build_setup(std::string pathtofile, bool verbose = true);
	void build_augmented(std::string pathtofile, bool verbose = true);
	virtual void build_extra(std::string pathtofile, bool verbose = true);
	void build(std::string params_file = UNDEFINED, std::string data_file = UNDEFINED,
			   std::string augm_file = UNDEFINED, std::string setup_file = UNDEFINED,
			   std::string extra_file = UNDEFINED, bool verbose = true);

	//===== Display Internal Information =====//
	void show_setup();
	void show_summary();
	void show_report();
	void show_parameters();
	void show_param_diagnostics();
	void show_results();

	//===== Return Internal Information with Rcpp =====//
#ifdef USE_RCPP
	Rcpp::List return_setup();
#endif

	//=====General Initialization=====//
	void set_rate_per_group(Individual &this_indiv);
	void set_rate_per_contact(int h);
	virtual void assign_rate_per_group_to_individual(Individual &this_indiv);
	virtual void assign_rate_per_contact_to_household(int h, int from, int to);
	virtual void initialize_model() {};
	virtual void initialize_chain() {};

	//=====Augmented Variables=====//
	void assign_augmented_variables();
	void init_augmented_variables();
	virtual void update_first_case(int h);
	virtual void augm_var_store_dependecies(int h, int i);
	virtual double augm_var_propose(int h, int i, std::string augmID);
	virtual void augm_var_update(int h, int i, std::string augmID, double new_val);
	virtual double get_infec_time(int h, int i) const;

	//=====Log-Likelihood computations=====//
	void set_global_loglike(long double x) { global_loglike = x; };
	long double get_global_loglike() const { return global_loglike; };
	virtual double compute_single_house_loglike(int h);
	virtual long double compute_global_loglike();
	long double run_eval_loglikelihood();

	//=====Monte Carlo Markov Chain moves=====//
	void init_move_counts();
	void single_param_move(std::string parID);
	void move_all_parameters();
	double single_household_augm_move(int h);
	double get_augm_accept_rate(std::string augmID);
	void move_all_augm_variables();
	void move_and_adapt_all_parameters(ulong step);
	void init_model_diagnostics();
	void update_model_diagnostics();

	//=====Infer=====//
	void run_inference(std::string output_filename, bool verbose = true, bool progbar = true);

	//=====Simulate=====//
	virtual void simulate_single_house(int h);
	virtual void format_augmented_data() {};
	void run_simulation(std::string output_filename, bool verbose = true);

	//=====Save Inference State=====//
	void init_inferstate_file(std::string output_filename);
	void print_inferstate();
	void close_inferstate_file() { inferstate_ofstream.close(); }

	//=====Save augmented variables=====//
	void print_augm_variables(std::string output_filename, bool verbose = true);

private:
	long double global_loglike = 0.0;

	//=====BUILD=====//

	//-----build : auxiliary attributes-----//
	std::ofstream inferstate_ofstream;			   // outoput inference state file
	std::vector<ssargs> load_augm_definition = {}; // temporary loader for augmented variable definition

	//----build : auxiliary functions-----//
	std::string string_cleaner(std::string token);
	bool stobool(std::string token);
	ssargs parse_dict_string(std::string token);
	sdargs parse_dict_args(std::string token);
	std::vector<std::string> read_column_headers(std::ifstream &input,
												 std::vector<std::string> required_headers = {},
												 std::string pathtofile = UNDEFINED);
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
