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
	double get_value(const std::string k) const { return params.at(k).get_value(); };
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); };
	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); };

	// Setter
	void set_value(const std::string k, double new_value) { params.at(k).set_value(new_value); };
	void init_value(const std::string k, std::mt19937_64 &rand_gen) { params.at(k).init_value(rand_gen); };
	void reset_moves(const std::string k) { params.at(k).reset_moves(); };
	void set_num_propos_pp(const std::string k) { params.at(k).set_num_propos_pp(); };
	void set_num_accept_pp(const std::string k) { params.at(k).set_num_accept_pp(); };
	void update_accept_rate(const std::string k, bool move_accepted, int step) { params.at(k).update_accept_rate(move_accepted, step); };
	void adapt(const std::string k, int step, int m_num_adapt_iter) { params.at(k).adapt(step, m_num_adapt_iter); };

	void add_param(const std::string k, Parameter new_param);

	// Computations
	bool valid_param(const std::string k, double new_value) { return params.at(k).valid(new_value); };
	double propose(const std::string k, std::mt19937_64 &rand_gen, double old_value) { return params.at(k).propose(rand_gen, old_value); };
	double log_ratio_proposal(const std::string k, double old_value, double new_value) { return params.at(k).propose_log_ratio(old_value, new_value); };
	double log_ratio_prior(const std::string k, double old_value, double new_value) { return params.at(k).prior_log_ratio(old_value, new_value); };

	// Display
	void print();
};

// global params declaration // FIXME : find a better position maybe
extern Parameter_Dictionary glob_params;

/////////////////////////////////////////////////////////
/////// Augmented : Class declaration /////////
////////////////////////////////////////////////////////

class Augmented : public Parameter
{
private:
	// Attributes
	std::map<std::string, std::string> model_arg_names = {}; // arguments to the model function

public:
	// Constructors
	// Constructors
	Augmented(
		double new_value,
		std::string strat_name = "_none_",
		kwargs strat_args = {},
		kwargs value_range = {},
		std::string model_name = "_none_",
		std::map<std::string, std::string> model_args = {}) : Parameter(new_value, strat_name, strat_args, value_range, model_name)
	{
		for (auto const &entry : model_args)
		{
			model_arg_names.emplace(entry.first, entry.second);
			args.emplace(entry.first, glob_params.get_value(entry.second));
		}
	};

	kwargs &address_args() override;
};

#endif // PARAMETER__H