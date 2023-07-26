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

#define MCMC_OPTIMAL_ACCEPTANCE 0.25
#define N_ITER_TIME_INFEC 1

//----------------------------------------------------------------------
// McmcObject : Class declaration
//----------------------------------------------------------------------

class McmcObject {
	// Attributes
	std::vector<std::string> m_param_names;
	std::vector<std::string> m_param_to_infer_names;
	std::mt19937_64 m_rand_gen;
	size_t m_iterations = 0;
	size_t m_iter_time_infec = 0;
	size_t m_adaptive_steps = 0;
	double m_global_LogLike = 0.0;
	std::vector<Household> m_data = {};
	std::vector<double> m_hh_LogLike = {};

	public:
		int m_num_accept_data = 0;
		int m_num_propos_data = 0;

		// Constructors
		McmcObject ();
		McmcObject (
			unsigned seed,
			size_t num_iterations,
			size_t num_adaptive_steps,
			std::vector<Household> data,
			size_t n_iter_time_infec = N_ITER_TIME_INFEC
			);
		~McmcObject () {};

		// Getter
		size_t num_households () const { return m_data.size(); };
		size_t num_iterations () const { return m_iterations; };
		size_t get_num_iter_time_infec () const {return m_iter_time_infec; };
		double global_LogLike () const { return m_global_LogLike; };
		double hh_LogLike (int index) const { return m_hh_LogLike[index]; };
		std::vector<double> hh_LogLike () const { return m_hh_LogLike; };

		// Setter
		void init_param_values ();
		void reset_moves ();
		void init_augm_data (bool display = false);
		void init_LogLike (bool display = false);

		// Computations
		void update_param (std::string parID, size_t step, bool display = false);
		void update_augm_infec_times (bool display = false);
		void update_augm_symptom_onset (bool display = false);

		void run (
			std::string output_file,
			size_t num_sweeps,
			bool display = false
			);
};

#endif // MCMC_OBJECT__H
