// metropolis.hpp //

#ifndef METROPOLIS__H
#define METROPOLIS__H

#include <string>
#include <vector>
#include <cmath>
#include <random>
#include <iostream>

#include "utils.hpp"

//----------------------------------------------------------------------
// Metropolis : class
//----------------------------------------------------------------------

// base class to deal with metropolis strategies
class Metropolis
{
	// strategies 
	double (Metropolis::*_propose_func_)(std::mt19937_64& rand_gen, double old_value);
	double (Metropolis::*_log_ratio_propose_func_)(double old_value, double new_value);
	void (Metropolis::*_adapt_func_)(int step, int num_adapt_steps);

	// valididation
	double _range_low_ = -INFINITY;
	double _range_up_ = INFINITY;
	bool (Metropolis::*_valid_func_)(double new_value);

public:
	// Attributes
	kwargs args = {};
	int num_accept = 0;
	int num_propos = 0;
	double accept_rate = 0.0;

	// Constructor
	Metropolis(
		std::string = "_none_",
		kwargs init_args = {},
		kwargs init_range = {});
	Metropolis() = default;
	~Metropolis(){};

	// Setter
	void reset_moves();
	void set_range(kwargs new_range);
	void set_args(kwargs new_args);
	void update_accept_rate(bool move_accepted, int step);

	// Computations
	bool valid(double new_value) { return (this->*_valid_func_)(new_value); };

	// Methods
	double propose(std::mt19937_64& rand_gen, double old_value) { return (this->*_propose_func_)(rand_gen, old_value); };
	double log_ratio_proposal(double old_value, double new_value) { return (this->*_log_ratio_propose_func_)(old_value, new_value); };
	void adapt(int step, int num_adapt_steps) { (this->*_adapt_func_)(step, num_adapt_steps); };

private:
	// --- Strategy : functional methods --- //

	// normal
	double _normal_propose_(std::mt19937_64& rand_gen, double old_value);
	double _normal_log_ratio_proposal_(double old_value, double new_value);
	void _normal_adapt_(int step, int num_adapt_iter);

	// exp
	double _exp_propose_(std::mt19937_64& rand_gen, double old_value);
	double _exp_log_ratio_proposal_(double old_value, double new_value);
	void _exp_adapt_(int step, int num_adapt_iter);

	// gamma
	double _gamma_propose_(std::mt19937_64& rand_gen, double old_value);
	double _gamma_log_ratio_proposal_(double old_value, double new_value);
	void _gamma_adapt_(int step, int num_adapt_iter);

	// --- Validation : functional methods --- //
	bool _valid_dummy_(double new_value);
	bool _valid_range_low_(double new_value);
	bool _valid_range_up_(double new_value);
	bool _valid_range_both_(double new_value);

};

#endif // METROPOLIS__H