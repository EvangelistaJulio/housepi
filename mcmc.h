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
// #include <omp.h>

#include "model.h"

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

public:
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> unif;
	Setup_MCMC_Struct setup = {};
	int m_num_accept_data = 0; // FIXME : to be cleaned out
	int m_num_propos_data = 0;

	// Constructors
	MCMC(Setup_MCMC_Struct setup) : unif(0.0, 1.0)
	{
		this->setup = setup;
		rng.seed(setup.seed);
	}
	~MCMC() {}

	// --- Getter --- //
	double get_global_loglike() const { return global_loglike; };

	// --- Setter --- //
	void reset_moves();
	void initialize_model(bool random_params = false);

	// --- Computations --- //
	void init_global_loglike();
	void move_parameters(double inv_temperature = 1.0);
	void move_augm_variables(double inv_temperature = 1.0);

	//=====RUN=====//
	void run(std::string output_filename, bool random_params = false, double inv_temperature = 1.0);

	//=====SIMULATE=====//
	void simulate(std::string output_filename);

	// --- Save parameters --- //
	void init_param_file(std::string output_filename);
	void print_parameters(unsigned step);
	void close_param_file() { param_ofstream.close(); }

	// --- Save augmented variables --- //
	void print_augm_variables(std::string output_filename);

private:
	// --- infer : attributes --- //
	std::ofstream param_ofstream;
};

#endif // MCMC__H
