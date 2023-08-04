// parameter.hpp //

#ifndef PARAMETER__H
#define PARAMETER__H

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <iostream>

#include "utils.hpp"

typedef double (*prior_pdf_ptr)(double x, std::vector<double> args);
typedef double (*prior_rand_ptr)(std::mt19937_64 &rand_gen, std::vector<double> args);

//----------------------------------------------------------------------
// Metropolis Strategy auxiliary classes
//----------------------------------------------------------------------

// base class to deal with metropolis strategies
class Metropolis_Strategy
{
	std::string name = "_none_";

public:
	// Attributes
	std::vector<double> range = {};
	double rate = 0.0;
	int num_accept = 0;
	int num_propos = 0;
	double accept = 0.0;

	// Constructor
	Metropolis_Strategy(double rate, std::vector<double> new_range);

	// Setter
	void reset_moves();
	void set_rate(double new_rate) { rate = new_rate; };
	void adapt_rate(bool move_accepted, int step, int num_adapt_iter);

	// Computations
	virtual double propose_move(std::mt19937_64 &rand_gen, double old_value) { return old_value + unif_rand(rand_gen); };
	virtual double log_ratio_proposal() { return 0.; };
	bool within_range(double new_value);
};

class Strategy_exp : public Metropolis_Strategy
{
	std::string name = "exp";

public:
	// Constructor
	Strategy_exp(double rate, std::vector<double> new_range) : Metropolis_Strategy(rate, new_range){};

	// Computations
	double propose_move(std::mt19937_64 &rand_gen, double new_value) override;
	double log_ratio_proposal() override;
};

class Strategy_normal : public Metropolis_Strategy
{
	std::string name = "normal";

public:
	// Constructor
	Strategy_normal(double rate, std::vector<double> new_range) : Metropolis_Strategy(rate, new_range){};

	// Computations
	double propose_move(std::mt19937_64 &rand_gen, double new_value) override;
	double log_ratio_proposal() override;
};

////////////////////////////////////////////////////////////
//////// Parameter : Class declaration /////////////////////
////////////////////////////////////////////////////////////

class Parameter
{
	// Attributes
	double p_value = 0.;					 // value of the parameter
	std::string p_prior = "_const_";		 // choice of the prior to the parameter
	std::vector<double> p_prior_args = {};	 // arguments to the prior function
	prior_pdf_ptr p_prior_log_pdf = nullptr; // pointer to the log prior pdf
	prior_rand_ptr p_prior_rand = nullptr;	 // the actual generator to execute

public:
	std::string p_name = "_no_name_"; // parameter name
	bool p_infer = false;			  // wheter to infer the parameter or not

	// Constructors
	Parameter(
		std::string name,
		double value,
		bool infer,
		std::string prior = "_const_",
		std::vector<double> args = {},
		std::string metropolis_strategy = "exp",
		double metropolis_rate = 0.,
		std::vector<double> range = {});
	Parameter() = default;
	~Parameter(){};

	// Getter
	double get_value() const;

	// Setter
	void set_value(double value);

	// Computations
	double generate_from_prior(std::mt19937_64 &rand_gen);
	double log_ratio_prior(double new_value);

	// --- Aliases for Metropolis Strategy ---
	// -> Getter
	int get_num_accept() const { return _metropolis->num_accept; }
	int get_num_propos() const { return _metropolis->num_propos; }
	double get_accept() const { return _metropolis->accept; }

	// -> Setter
	void set_rate_metropolis(double rate);
	void reset_moves() { _metropolis->reset_moves(); }
	void set_num_accept_pp() { _metropolis->num_accept++; }
	void set_num_propos_pp() { _metropolis->num_propos++; }
	void set_accept(double accept) { _metropolis->accept = accept; }
	void adapt_rate(bool move_accepted, int step, int num_adapt_iter) { _metropolis->adapt_rate(move_accepted, step, num_adapt_iter); }

	// -> Computations
	double propose_move(std::mt19937_64 &rand_gen) { return _metropolis->propose_move(rand_gen, p_value); }
	bool within_range(double new_value) { return _metropolis->within_range(new_value); }
	double log_ratio_proposal() { return _metropolis->log_ratio_proposal(); }

private:
	Metropolis_Strategy *_metropolis; // Pointer for polymorphism
};

///////////////////////////////////////////////////////////
//////// Par_Dict : Class declaration /////////////////////
///////////////////////////////////////////////////////////

// alias class for global parameters
class Par_Dict
{
	// attributes
	std::map<std::string, Parameter> parameters = {};

public:
	std::vector<std::string> to_be_inferred = {};

	// Getter
	double value(std::string parID) const;
	bool infer(std::string parID) const;
	int num_propos(std::string parID) const;
	int num_accept(std::string parID) const;

	// Setter
	void add_param(std::string parID, Parameter new_param);
	void set_value(std::string parID, double new_value);
	void init_value(std::string parID, std::mt19937_64 &m_rand_gen);
	void reset_moves(std::string parID);
	void set_num_propos_pp(std::string parID);
	void set_num_accept_pp(std::string parID);
	void adapt_rate(std::string parID, bool move_accepted, int step, int m_num_adapt_iter);

	// Computations
	bool within_range(std::string parID, double new_value);
	double propose_move(std::string parID, std::mt19937_64 &m_rand_gen);
	double log_ratio_proposal(std::string parID);
	double log_ratio_prior(std::string parID, double new_value);

	// Other
	void print();
};

// global parameters declaration
extern Par_Dict glob_params;

#endif // PARAMETER__H