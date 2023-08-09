// parameter.hpp //

#ifndef PARAMETER__H
#define PARAMETER__H

#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <iostream>

#include "metropolis.hpp"

typedef double (*pdf_ptr)(double x, kwargs args);
typedef double (*rand_ptr)(std::mt19937_64 &rand_gen, kwargs args);

////////////////////////////////////////////////////////////
//////// Parameter : Class declaration /////////////////////
////////////////////////////////////////////////////////////

// FIXME : this should be inherited from metropolis :
// most of the methods comes from it; it should have value as attribute
// in order to be expanded to augmented variables also

class Parameter
{
	// Attributes
	double p_value = 0.;					 	// value of the parameter
	std::string p_prior = "_const_";			 // choice of the prior to the parameter
	kwargs p_prior_args = {};	 				// arguments to the prior function
	pdf_ptr p_prior_log_pdf = nullptr; 	// pointer to the log prior pdf
	rand_ptr p_prior_rand = nullptr;	 	// the actual generator to execute

public:
	std::string p_name = "_none_"; 			// parameter name
	bool p_infer = false;			 		 // wheter to infer the parameter or not

	// Constructors
	Parameter(
		std::string name,
		double value,
		bool infer,
		std::string prior = "_const_",
		kwargs args = {},
		std::string metropolis_strategy = "exp",
		double metropolis_rate = 0.,
		std::vector<double> range = {});
	Parameter() = default;
	~Parameter(){};

	// Getter
	double get_value() const { return p_value; };

	// Setter
	void set_value(double new_value) { p_value = new_value; };

	// Computations
	double generate_from_prior(std::mt19937_64 &rand_gen);
	double log_ratio_prior(double new_value);

	// --- Aliases for Metropolis Strategy ---
	// -> Getter
	int get_num_accept() const { return _tmp_->num_accept; }
	int get_num_propos() const { return _tmp_->num_propos; }
	double get_accept() const { return _tmp_->accept; }

	// -> Setter
	void set_rate_tmp_(double new_rate) { _tmp_->set_args({{"rate", new_rate}}); };
	void reset_moves() { _tmp_->reset_moves(); }
	void set_num_accept_pp() { _tmp_->num_accept++; }
	void set_num_propos_pp() { _tmp_->num_propos++; }
	void set_accept(double accept) { _tmp_->accept = accept; }
	void adapt_rate(bool move_accepted, int step, int num_adapt_iter) { _tmp_->adapt_rate(move_accepted, step, num_adapt_iter); }

	// -> Computations
	double propose(std::mt19937_64 &rand_gen) { return _tmp_->propose(rand_gen, p_value); }
	bool within_range(double new_value) { return _tmp_->within_range(new_value); }
	double log_ratio_proposal() { return _tmp_->log_ratio_proposal(); }

private:
	Strategy_Param *_tmp_; // Pointer for polymorphism
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
	double propose(std::string parID, std::mt19937_64 &m_rand_gen);
	double log_ratio_proposal(std::string parID);
	double log_ratio_prior(std::string parID, double new_value);

	// Other
	void print();
};

// global parameters declaration
extern Par_Dict glob_params;

#endif // PARAMETER__H