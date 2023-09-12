// mcmc.hpp //

#ifndef MCMC_OBJECT__H
#define MCMC_OBJECT__H

#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
//#include <omp.h>

#include "household.hpp"

// Auxiliary struct to load all setup information
struct Setup_Struct
{
public:
	// Attributes
	unsigned seed = 0;
	unsigned num_iter = 0;
	unsigned num_adapt_iter = 0;
	unsigned num_augm_sweeps = 0;
	unsigned thinning = 0;
	char delim = ' ';
	char dict_delim = '/';
	bool verbose = false;

	// FIXME : attributes toi be updated
	bool is_data_augm = false;
};

//----------------------------------------------------------------------
// MCMC_Class : Class declaration
//----------------------------------------------------------------------

class MCMC_Class
{
private:
	std::mt19937_64 m_rng;
	double m_global_loglike = 0.0;
	std::vector<Household> m_data = {};
	std::uniform_real_distribution<double> m_unif;
	Setup_Struct m_setup = {};

public:
	int m_num_accept_data = 0;
	int m_num_propos_data = 0;
	Onset_to_Transm infectivity_profile;

	// Constructors
	MCMC_Class(
		Setup_Struct setup,
		std::vector<Household> data,
		kwargs_ext infec_prof_args) : m_unif(0.0, 1.0), infectivity_profile(infec_prof_args)
	{
		m_setup = setup;
		m_rng.seed(setup.seed);
		m_data = data;
	}
	~MCMC_Class(){};

	// --- Getter --- //
	size_t num_households() const { return m_data.size(); };
	double get_global_loglike() const { return m_global_loglike; };

	// --- Setter --- //
	void reset_moves();
	void init_augm_data(std::mt19937_64 &rng);

	// --- Computations --- //
	void update_global_loglike();
	void update_parameters();
	void update_augm_variables();

	// --- INFER --- //
	void infer(std::string output_filename);

	// --- Save parameters --- //
	void init_param_file(std::string output_filename);
	void print_parameters(unsigned step);
	void close_param_file() { param_ofstream.close(); }

	// --- Save augmented variables --- //
	void open_augm_file(std::string output_filename);
	void print_augm_variables();
	void close_augm_file() { augm_ofstream.close(); }

private:
	// --- infer : attributes --- //
	std::ofstream param_ofstream;
	std::ofstream augm_ofstream;

	void progress_bar(size_t iter, size_t n_iter);
};

#endif // MCMC_OBJECT__H
