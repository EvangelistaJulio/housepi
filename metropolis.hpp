// metropolis.hpp //

#ifndef METROPOLIS__H
#define METROPOLIS__H

#include <string>
#include <vector>
#include <map>
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
	std::string name = "_none_";

public:
	// Attributes
	std::map<std::string, double> args = {};
	std::vector<double> range = {};
	int num_accept = 0;
	int num_propos = 0;
	double accept = 0.0;

	// Constructor
	Metropolis(std::map<std::string, double> init_args, std::vector<double> init_range);

	// Setter
	void reset_moves();
	void set_range(std::vector<double> new_range);
	void set_args(std::map<std::string, double> new_args);

	// Computations
	bool within_range(double new_value);

	// Virtual methods
	virtual double propose(std::mt19937_64 &rand_gen, double old_value) { return old_value + unif_rand(rand_gen); };
	virtual double log_ratio_proposal() { return 0.; };
};

//----------------------------------------------------------------------
// Metropolis : Parameter Strategy auxiliary classes
//----------------------------------------------------------------------

class Strategy_Param : public Metropolis
{
	std::string name = "_none_";

public:
	// Constructor
	Strategy_Param(double init_rate, std::vector<double> init_range) : Metropolis({{"rate", init_rate}}, init_range){};

	void adapt_rate(bool move_accepted, int step, int num_adapt_iter);
};

class Strategy_Param_exp : public Strategy_Param
{
	std::string name = "exp";

public:
	// Constructor
	Strategy_Param_exp(double init_rate, std::vector<double> init_range) : Strategy_Param(init_rate, init_range){};

	// Computations
	double propose(std::mt19937_64 &rand_gen, double new_value) override;
	double log_ratio_proposal() override;
};

// normal metropolis stratregy
class Strategy_Param_normal : public Strategy_Param
{
	std::string name = "normal";

public:
	// Constructor
	Strategy_Param_normal(double init_rate, std::vector<double> init_range) : Strategy_Param(init_rate, init_range){};

	// Computations
	double propose(std::mt19937_64 &rand_gen, double new_value) override;
	double log_ratio_proposal() override;
};

#endif // METROPOLIS__H