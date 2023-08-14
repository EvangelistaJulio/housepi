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

class Parameter : public Metropolis
{
private:
	// Attributes
	kwargs p_prior_args = {}; // arguments to the prior function

public:
	// FIXME : prior_args !!!

	// Constructors
	using Metropolis::Metropolis;
};

///////////////////////////////////////////////////////////
/////// Parameter_Dictionary : class declaration //////////
///////////////////////////////////////////////////////////

// This is a wrapper class containing aliases for the global params.
class Parameter_Dictionary
{
private:
	// Attributes
	std::map<std::string, Parameter> params = {};

public:
	// Attributes
	std::vector<std::string> to_be_inferred = {};

	// Getter
	double get_value(std::string k) const { return params.at(k).get_value(); };
	bool get_infer(std::string k) const { return params.at(k).get_infer(); };
	double get_accept_rate(std::string k) const { return params.at(k).get_accept_rate(); };

	// Setter
	void set_value(std::string k, double new_value) { params.at(k).set_value(new_value); };
	void init_value(std::string k, std::mt19937_64 &rand_gen) { params.at(k).init_value(rand_gen); };
	void reset_moves(std::string k) { params.at(k).reset_moves(); };
	void set_num_propos_pp(std::string k) { params.at(k).set_num_propos_pp(); };
	void set_num_accept_pp(std::string k) { params.at(k).set_num_accept_pp(); };
	void update_accept_rate(std::string k, bool move_accepted, int step) { params.at(k).update_accept_rate(move_accepted, step); };
	void adapt(std::string k, int step, int m_num_adapt_iter) { params.at(k).adapt(step, m_num_adapt_iter); };

	void add_param(std::string k, Parameter new_param);

	// Computations
	bool valid_param(std::string k, double new_value) { return params.at(k).valid(new_value); };
	double propose(std::string k, std::mt19937_64 &rand_gen, double old_value) { return params.at(k).propose(rand_gen, old_value); };
	double log_ratio_proposal(std::string k, double old_value, double new_value) { return params.at(k).propose_log_ratio(old_value, new_value); };
	double log_ratio_prior(std::string k, double old_value, double new_value) { return params.at(k).prior_log_ratio(old_value, new_value); };

	// Display
	void print();
};


// global params declaration // FIXME : find a better position maybe
extern Parameter_Dictionary glob_params;

/////////////////////////////////////////////////////////
/////// Augmented : Class declaration /////////
////////////////////////////////////////////////////////

class Augmented : public Metropolis
{
private:
	// Attributes
	std::map<std::string, std::string> a_prior_args = {}; // arguments to the prior function

public:
	// Constructors
	using Metropolis::Metropolis;

	// Computations
	kwargs get_prior_args() override;
};

#endif // PARAMETER__H