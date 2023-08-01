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

#define NUM_AUGM_ITER 1
#define DEFAULT_SEED 20210415

//----------------------------------------------------------------------
// MCMC_Class : Class declaration
//----------------------------------------------------------------------

class MCMC_Class {
	// Attributes
	std::vector<std::string> m_param_to_infer_names;
	std::mt19937_64 m_rand_gen;
	int m_num_iter = 0;
	int m_num_augm_iter = 0;
	int m_num_adapt_iter = 0;
	double m_global_LogLike = 0.0;
	std::vector<Household> m_data = {};
	std::vector<double> m_hh_LogLike = {};

	public:
		int m_num_accept_data = 0;
		int m_num_propos_data = 0;

		// Constructors
		MCMC_Class ();
		MCMC_Class (
			unsigned seed,
			int num_num_iter,
			int num_adapt_iter,
			std::vector<Household> data,
			int num_augm_iter = NUM_AUGM_ITER
			);
		~MCMC_Class () {};

		// Getter
		size_t num_households () const { return m_data.size(); };
		int get_num_iter () const { return m_num_iter; };
		int get_num_augm_iter () const {return m_num_augm_iter; };
		double global_LogLike () const { return m_global_LogLike; };
		double hh_LogLike (int index) const { return m_hh_LogLike[index]; };
		std::vector<double> hh_LogLike () const { return m_hh_LogLike; };

		// Setter
		void init_param_values ();
		void reset_moves ();
		void init_augm_data ();
		void init_LogLike ();

		// Computations
		bool update_param (std::string parID);
		void update_augm_infec_times ();
		void update_augm_symptom_onset ();

		void run (
			std::string output_filename,
			int num_sweeps
			);
};

#endif // MCMC_OBJECT__H
