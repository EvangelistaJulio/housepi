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

#include "household.hpp"

//----------------------------------------------------------------------
// MCMC_Class : Class declaration
//----------------------------------------------------------------------

class MCMC_Class
{
private:
	std::mt19937_64 m_rng;
	double m_global_loglike = 0.0;
	std::vector<Household> m_data = {};
	std::uniform_real_distribution<double> unif;

public:
	int m_num_accept_data = 0;
	int m_num_propos_data = 0;
	Onset_to_Transm infectivity_profile;

	// Constructors
	MCMC_Class(
		unsigned seed,
		std::vector<Household> data,
		kwargs_ext infec_prof_args) : infectivity_profile(infec_prof_args)
	{
		m_rng.seed(seed);
		m_data = data;
	}
	~MCMC_Class(){};

	// Getter
	size_t num_households() const { return m_data.size(); };
	double get_global_loglike() const { return m_global_loglike; };

	// Setter
	void reset_moves();
	void init_augm_data(std::mt19937_64 &rng);
	void init_global_loglike();

	// Computations
	bool update_param(std::string parID);
	void update_augm_variables();

	void infer(
		std::string output_filename,
		int num_iter = DEFAULT_NUM_ITER,
		int num_param_sweeps = DEFAULT_NUM_PARAM_SWEEEPS,
		int num_augm_sweeps = DEFAULT_NUM_AUGM_SWEEPS,
		int num_adapt_steps = DEFAULT_NUM_ADAPT_STEPS,
		std::string tab = TAB);
};

#endif // MCMC_OBJECT__H
