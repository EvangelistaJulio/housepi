// mcmc.hpp //

#ifndef MCMC_OBJECT__H
#define MCMC_OBJECT__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>

#include "household.hpp"

//----------------------------------------------------------------------
// MCMC_Class : Class declaration
//----------------------------------------------------------------------

class MCMC_Class
{
	// Attributes
	std::mt19937_64 m_rand_gen;
	double m_global_loglike = 0.0;
	std::vector<Household> m_data = {};

public:
	int m_num_accept_data = 0;
	int m_num_propos_data = 0;

	// Constructors
	MCMC_Class();
	MCMC_Class(
		unsigned seed,
		std::vector<Household> data);
	~MCMC_Class(){};

	// Getter
	size_t num_households() const { return m_data.size(); };
	double get_global_loglike() const { return m_global_loglike; };

	// Setter
	void init_param_to_infer();
	void reset_moves();
	void init_augm_data();
	void init_global_loglike();

	// Computations
	bool update_param(std::string parID);
	void update_augm_variables();

	void infer(
		std::string output_filename,
		int num_iter = DEFAULT_NUM_ITER,
		int num_param_sweeps = DEFAULT_NUM_PARAM_SWEEEPS,
		int num_augm_sweeps = DEFAULT_NUM_AUGM_SWEEPS,
		int num_adapt_iter = DEFAULT_NUM_ADAPT,
		std::string tab = TAB);
};

#endif // MCMC_OBJECT__H
