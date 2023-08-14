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
private:
	// general attributes
	double value;			// value of the parameter
	double value_0;			// initial value
	bool infer = false; 	// inference mode

	// prior
	pdf_ptr _prior_pdf_;	 // pointer to the prior pdf
	pdf_ptr _prior_log_pdf_; // pointer to the log prior pdf
	rand_ptr _prior_rand_;	 // the actual generator to execute

	// strategies
	kwargs propose_args = {};
	double (Metropolis::*_propose_func_)(std::mt19937_64 &rand_gen, double old_value);
	double (Metropolis::*_propose_log_ratio_func_)(double old_value, double new_value);
	void (Metropolis::*_adapt_func_)(int step, int num_adapt_steps);

	// valididation
	double _range_low_ = -INFINITY;
	double _range_up_ = INFINITY;
	bool (Metropolis::*_valid_func_)(double new_value);

	// acceptation
	int num_accept = 0;
	int num_propos = 0;
	double accept_rate = 0.0;

public:

	// Constructor
	Metropolis(
		double new_value,
		std::string strat_name = "_none_",
		kwargs strat_args = {},
		kwargs value_range = {},
		std::string prior_name = "_none_",
		kwargs prior_args = {});
	Metropolis() = default;
	~Metropolis(){};

	// Getter
	double get_value() const { return value; };
	double get_value_0() const { return value_0; };
	bool get_infer() const { return infer; };
	double get_accept_rate() const { return accept_rate; };
	virtual kwargs get_prior_args();

	// Setter
	void set_value(double new_value) { value = new_value; };
	void init_value(std::mt19937_64 &rand_gen) { value = prior_generate(rand_gen); };
	void set_range(kwargs new_range);
	void set_propose_args(kwargs new_args);
	void update_accept_rate(bool move_accepted, int step);
	void set_num_accept_pp() { num_accept++; };
	void set_num_propos_pp() { num_propos++; };
	void reset_moves();

	// Computations
	bool valid(double new_value) { return (this->*_valid_func_)(new_value); };

	// Methods for strategies
	double propose(std::mt19937_64 &rand_gen, double old_value) { return (this->*_propose_func_)(rand_gen, old_value); };
	double propose_log_ratio(double old_value, double new_value) { return (this->*_propose_log_ratio_func_)(old_value, new_value); };
	void adapt(int step, int num_adapt_steps) { (this->*_adapt_func_)(step, num_adapt_steps); };

	// Methods for prior
	double prior_generate(std::mt19937_64 &rand_gen) { return _prior_propose_(rand_gen, 0.); }
	double prior_log_ratio(double old_value, double new_value) { return _prior_propose_log_ratio_(old_value, new_value); }

private:
	// --- Prior : functional methods --- //

	double _prior_propose_(std::mt19937_64 &rand_gen, double _);
	double _prior_propose_log_ratio_(double old_value, double new_value);

	// --- Strategy : functional methods --- //

	// normal
	double _normal_propose_(std::mt19937_64 &rand_gen, double old_value);
	double _normal_propose_log_ratio_(double old_value, double new_value);
	void _normal_adapt_(int step, int num_adapt_iter);

	// lognorm
	double _lognorm_propose_(std::mt19937_64 &rand_gen, double old_value);
	double _lognorm_propose_log_ratio_(double old_value, double new_value);
	void _lognorm_adapt_(int step, int num_adapt_iter);

	// gamma
	double _gamma_propose_(std::mt19937_64 &rand_gen, double old_value);
	double _gamma_propose_log_ratio_(double old_value, double new_value);
	void _gamma_adapt_(int step, int num_adapt_iter);

	// exp
	double _exp_propose_(std::mt19937_64 &rand_gen, double old_value);
	double _exp_propose_log_ratio_(double old_value, double new_value);
	void _exp_adapt_(int step, int num_adapt_iter);

	// --- MCMC moves --- //
	void _adapter_(std::string argID, int step, int num_adapt_iter);

	// --- Validation : functional methods --- //
	bool _valid_dummy_(double new_value);
	bool _valid_range_low_(double new_value);
	bool _valid_range_up_(double new_value);
	bool _valid_range_both_(double new_value);
};

#endif // METROPOLIS__H