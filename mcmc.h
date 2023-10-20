// mcmc.h //

#ifndef MCMC__H
#define MCMC__H

#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>

#include "model.h"

extern bool Verbose;

// Auxiliary struct to load all setup information
struct Setup_MCMC_Struct
{
public:
	// Attributes
	unsigned seed = 20210415;
	unsigned num_iter = 20000;
	unsigned num_adapt_iter = 5000;
	unsigned thinning = 20;
	timetype delta_t = 0.01;
	char delim = ' ';
	char dict_delim = '/';
};

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : Class ////////////////////////// //
////////////////////////////////////////////////////////////////////////

class MCMC : public Model
{
private:
	double global_loglike = 0.0;
	std::ofstream param_ofstream;
	
public:
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> unif;
	Setup_MCMC_Struct setup = {};
	
	std::unordered_map<std::string, unsigned> m_num_accept_data = {};
	std::unordered_map<std::string, unsigned> m_num_propos_data = {};

	// Constructors
	MCMC(Setup_MCMC_Struct new_setup) : unif(0.0, 1.0)
	{
		setup = new_setup;
		rng.seed(setup.seed);
	}
	~MCMC() {}

	//=====Getter=====//
	double get_global_loglike() const { return global_loglike; };

	//=====Setter=====//
	void reset_moves();
	void set_global_loglike(double x) { global_loglike = x; };
	void initialize_model(bool is_rand_params=false);

	//=====Open data=====//
	void build_data(std::string data_pathtofile);
	void build_parameters(std::string data_pathtofile);
	void add_individual(Individual this_indiv);
	int assign_community_structure(std::string HH_index_column);

	//=====Computations=====//
	double compute_global_loglike();
	void move_parameters(double inv_temperature = 1.0);
	void move_augm_variables(double inv_temperature = 1.0);

	//=====INFER=====//
	void inference(std::string output_filename, bool is_rand_params = false, double inv_temperature = 1.0);

	//=====SIMULATE=====//
	void simulation(std::string output_filename);

	//=====Save parameters=====//
	void init_param_file(std::string output_filename);
	void print_parameters(unsigned step);
	void close_param_file() { param_ofstream.close(); }

	//=====Save augmented variables=====//
	void print_augm_variables(std::string output_filename);
};

#endif // MCMC__H
