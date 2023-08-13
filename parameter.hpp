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

////////////////////////////////////////////////////////////
//////// Parameter : class declaration /////////////////////
////////////////////////////////////////////////////////////

class Parameter
{
private:
	pdf_ptr p_prior_log_pdf = nullptr; // pointer to the log prior pdf
	rand_ptr p_prior_rand = nullptr;   // the actual generator to execute

public:
	// Attributes
	double p_value;			  // value of the parameter
	kwargs p_prior_args = {}; // arguments to the prior function
	bool p_infer = false;	// wheter to infer the parameter or not
	Metropolis *metropolis; // Pointer for polymorphism

	// Constructors
	Parameter(
		double value,
		bool infer,
		std::string metropolis_strategy = "_none_",
		kwargs metropolis_args = {},
		kwargs range = {},
		std::string prior_name = "_none_",
		kwargs prior_args = {});
	Parameter() = default;
	~Parameter(){};

	// Getter
	double get_value() const { return p_value; };
	void set_value(double new_value) { p_value = new_value; };

	// Computations
	double generate_from_prior(std::mt19937_64& rand_gen);
	double log_ratio_prior(double old_value, double new_value);

	// --- Aliases from Metropolis Strategy ---
	// -> Getter
	int get_num_accept() const { return metropolis->num_accept; }
	int get_num_propos() const { return metropolis->num_propos; }
	double get_accept_rate() const { return metropolis->accept_rate; }

	// -> Setter
	void set_args(kwargs new_args) { metropolis->set_args(new_args); };
	void reset_moves() { metropolis->reset_moves(); }
	void set_num_accept_pp() { metropolis->num_accept++; }
	void set_num_propos_pp() { metropolis->num_propos++; }
	void update_accept_rate(bool move_accepted, int step) { metropolis->update_accept_rate(move_accepted, step); }
	void adapt(int step, int num_adapt_iter) { metropolis->adapt(step, num_adapt_iter); }

	// -> Computations
	bool valid(double new_value) { return metropolis->valid(new_value); }
	double log_ratio_proposal(double old_value, double new_value) { return metropolis->log_ratio_proposal(old_value, new_value); }
	double propose(std::mt19937_64& rand_gen, double old_value) { return metropolis->propose(rand_gen, old_value); }
};

///////////////////////////////////////////////////////////
/////// Parameter_Dictionary : class declaration //////////
///////////////////////////////////////////////////////////

// This is a wrapper class containing aliases for the global parameters.
class Parameter_Dictionary
{
private:
	// Attributes
	std::map<std::string, Parameter> parameters = {};

public:
	// Attributes
	std::vector<std::string> to_be_inferred = {};

	// Getter
	double get_value(std::string parID) const;
	bool get_infer(std::string parID) const;
	double get_accept_rate(std::string parID) const;

	// Setter
	void add_param(std::string parID, Parameter new_param);
	void set_value(std::string parID, double new_value);
	void init_value(std::string parID, std::mt19937_64 &m_rand_gen);
	void reset_moves(std::string parID);
	void set_num_propos_pp(std::string parID);
	void set_num_accept_pp(std::string parID);
	void update_accept_rate(std::string parID, bool move_accepted, int step);
	void adapt(std::string parID, int step, int m_num_adapt_iter);

	// Computations
	bool valid_param(std::string parID, double new_value);
	double propose(std::string parID, std::mt19937_64 &m_rand_gen, double old_value);
	double log_ratio_proposal(std::string parID, double old_value, double new_value);
	double log_ratio_prior(std::string parID, double old_value, double new_value);

	// Display
	void print();
};

// global parameters declaration
extern Parameter_Dictionary glob_params;

#endif // PARAMETER__H