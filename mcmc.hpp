// mcmc.hpp //

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
//#include <omp.h>

#include "model.hpp"

// Auxiliary struct to load all setup information
struct Setup_MCMC_Struct
{
public:
	// Attributes
	unsigned seed = 0;
	unsigned num_iter = 0;
	unsigned num_adapt_iter = 0;
	unsigned thinning = 0;
	bool is_data_augm = false;
	bool verbose = false;
	char delim = ' ';
	char dict_delim = '/';
};

////////////////////////////////////////////////////////////////////////
// ////////////////////////// MCMC : Class ////////////////////////// //
////////////////////////////////////////////////////////////////////////

class MCMC : public Model
{
private:
	std::mt19937_64 m_rng;
	double m_global_loglike = 0.0;
	std::vector<Household> m_data = {};
	std::uniform_real_distribution<double> m_unif;
	Setup_MCMC_Struct m_setup = {};

public:
	int m_num_accept_data = 0;
	int m_num_propos_data = 0;

	// Constructors
	MCMC(
		Setup_MCMC_Struct setup,
		std::vector<Household> data) : m_unif(0.0, 1.0)
	{
		m_setup = setup;
		m_rng.seed(setup.seed);
		m_data = data;
	}
	~MCMC(){}

	// --- Getter --- //
	size_t num_households() const { return m_data.size(); };
	double get_global_loglike() const { return m_global_loglike; };

	// --- Setter --- //
	void reset_moves();
	void declare_augm_data(std::mt19937_64 &rng);

	// --- Computations --- //
	void update_global_loglike();
	void move_parameters(double inv_temperature = 1.0);
	void move_augm_variables(double inv_temperature = 1.0);

	// --- INFER --- //
	void run(std::string output_filename, double inv_temperature = 1.0);

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
};

#endif // MCMC__H
