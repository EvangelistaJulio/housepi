#ifndef MCMCOBJECT__H
#define MCMCOBJECT__H

#include <random>
#include <vector>
#include "household.hpp"

// Class
class McmcObject {
	public:
		// Constructors
		McmcObject();
		McmcObject(
			unsigned seed,
			int num_iterations,
			int num_adaptive_steps,
			std::vector<Household> data,
			std::vector<double> parameter,
			std::vector<int> selected_parameter,
			std::vector<double> rate_for_random_walk,
			std::vector<int> lnormPrior,
			std::vector<double> sdLNormPrior,
			int n_iter_time_infec = 1, 
			double max_PCR_detectability = 10.0
			);
		~McmcObject() {};

		int iteration() const { return m_iterations; };
		int get_num_iter_time_infec() const {return m_iter_time_infec; };
		int proposed_move(int index) const {return m_num_proposed_move[index]; };
		int acceptedMove(int index) const {return m_num_accepted_move[index]; };
		int proposed_move_data() const {return m_proposed_move_data; };
		int accepted_move_data() const {return m_accepted_move_data; };
		double globalLogLik() const { return m_globalLogLik; };
		double hhLogLik(int index) const { return m_hhLogLik[index]; };
		std::vector<double> hhLogLik() const { return m_hhLogLik; };
		double rateRandomWalk(int index) const { return m_rate_for_random_walk[index]; };
		double parameter(int index) const { return m_parameter[index]; };
		size_t num_params() const { return m_parameter.size(); };
		
		void reset_moves();
		void init_param_values();
		void init_augm_data(bool display = false);
		void init_log_lik(bool display = false);
		void update_parameter(int parID, double step, bool display);
		void update_augm_infec_times(bool display);
		void update_augm_symptom_onset(bool display);

	private:
		// Attributes
		int m_iterations;
		int m_iter_time_infec;
		double m_adaptive_steps;
		std::mt19937_64 m_rnd_generator;
		double m_globalLogLik;
		std::vector<double> m_hhLogLik;
		int m_accepted_move_data;
		int m_proposed_move_data;
		std::vector<int> m_num_accepted_move;
		std::vector<int> m_num_proposed_move;
		std::vector<double> m_accept;
		std::vector<double> m_parameter;
		std::vector<int>  m_selected_param;
		std::vector<double> m_rate_for_random_walk;
		std::vector<int> m_lnormPrior;
		std::vector<double> m_sdLNormPrior;
		double m_max_PCR_detectability;
		double const m_opt_acceptance = 0.25;
		std::vector<Household> m_data;
		size_t m_nHH;
};

#endif
