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

///////////////////////////////////////////////////////////////////////
// ////////////////// Metropolis_Hastings : class ////////////////// //
///////////////////////////////////////////////////////////////////////

// FIXME : certain argument must be equal to 0.

class Metropolis_Hastings
{
protected:
	// acceptation
	int num_accept = 0;
	int num_propos = 0;
	double accept_rate = 0.0;

	double _adapter_(int step, int num_adapt_iter);

public:
	// Constructors
	Metropolis_Hastings() = default;
	~Metropolis_Hastings() {}

	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted, int step);
	void set_num_accept_pp() { num_accept++; };
	void set_num_propos_pp() { num_propos++; };
	void reset_moves();

	void set_range(kwargs new_range);
	bool valid(double new_value) { return (this->*_valid_func_)(new_value); }

	// --- functions (virtual) --- //
	virtual double propose(std::mt19937_64 &rng, double old_value)
	{
		std::uniform_real_distribution<double> distr_gen(-0.5, 0.5);
		return old_value + distr_gen(rng);
	}
	virtual double log_ratio() { return 0.0; } // NOTE : the options at the moment do not depend on old/new values
	virtual void adapt(int step, int num_adapt_iter)
	{
		(void)step;
		(void)num_adapt_iter;
	};

private:
	// --- Validation : functional methods --- //
	double _range_low_ = -INFINITY;
	double _range_up_ = INFINITY;
	bool (Metropolis_Hastings::*_valid_func_)(double new_value);

	bool _valid_dummy_(double new_value)
	{
		(void)new_value;
		return true;
	}
	bool _valid_range_low_(double new_value) { return (new_value > _range_low_); }
	bool _valid_range_up_(double new_value) { return (new_value < _range_up_); }
	bool _valid_range_both_(double new_value) { return _valid_range_low_(new_value) && _valid_range_up_(new_value); }
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Normal_Proposal : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Normal_Proposal : public Metropolis_Hastings
{
private:
	Normal_Distr distrib;
	std::vector<std::string> args_to_adapt = {"sd"};

public:
	// Constructors
	Normal_Proposal(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : distrib(new_args){};
	~Normal_Proposal() {}

	// --- functions (override) --- //
	double propose(std::mt19937_64 &rng, double old_value) override { return old_value + distrib.rand(rng); }
	double log_ratio() override { return 0.0; }
	void adapt(int step, int num_adapt_iter) override
	{
		for (std::string argID : args_to_adapt)
		{
			distrib.set_arg(argID, distrib.get_arg(argID) * _adapter_(step, num_adapt_iter));
		}
	}
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Lognorm_Proposal : class ////////////////// //
///////////////////////////////////////////////////////////////////////

class Lognorm_Proposal : public Metropolis_Hastings
{
private:
	Lognorm_Distr distrib;
	std::vector<std::string> args_to_adapt = {"sd"};

public:
	// Constructors
	Lognorm_Proposal(kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}}) : distrib(new_args){};
	~Lognorm_Proposal() {}

	// --- functions (override) --- //
	double propose(std::mt19937_64 &rng, double old_value) override { return old_value * distrib.rand(rng); }
	double log_ratio() override { return distrib.get_arg("sd"); }
	void adapt(int step, int num_adapt_iter) override
	{
		for (std::string argID : args_to_adapt)
		{
			distrib.set_arg(argID, distrib.get_arg(argID) * _adapter_(step, num_adapt_iter));
		}
	}
};

///////////////////////////////////////////////////////////////////////
// ////////////////////// Prior_Proposal : class /////////////////// //
///////////////////////////////////////////////////////////////////////

// FIXME : missing...

/////////////////////////////////////////////////////////////////
// ///////////////////// Parameter : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Parameter
{
private:
	// general attributes
	double *value_ptr;	// value of the parameter
	double value_0;		// initial value
	bool infer = false; // inference mode
	
protected:
	Distribution_1D _distribution_;
	Metropolis_Hastings _proposer_;

public:
	// Constructor
	Parameter(
		double *new_value_ptr,
		std::string strat_name = "_none_",
		kwargs strat_args = {},
		kwargs value_range = {},
		std::string prior_name = "_none_",
		kwargs prior_args = {});
	~Parameter() {}

	// Getter
	double get_value() const { return *(value_ptr); }
	double *get_address() const { return value_ptr; } // FIXME
	double get_value_0() const { return value_0; }
	bool get_infer() const { return infer; }

	// Setter
	void set_value(double new_value) { *(value_ptr) = new_value; }
	void init_value(std::mt19937_64 &rng) { *(value_ptr) = _distribution_.rand(rng); }

	// Computations
	bool valid(double new_value) { return _proposer_.valid(new_value); }

	// Methods for prior
	double prior_log_ratio(double old_value, double new_value) { return _distribution_.log_ratio(new_value, old_value); }

	// Methods for strategies
	double propose(std::mt19937_64 &rng, double old_value) { return (_proposer_.propose)(rng, old_value); }
	double propose_log_ratio() { return (_proposer_.log_ratio)(); }

	double get_accept_rate() const { return _proposer_.get_accept_rate(); }
	void set_num_propos_pp() { _proposer_.set_num_propos_pp(); }
	void set_num_accept_pp() { _proposer_.set_num_accept_pp(); }
	void update_accept_rate(bool is_accepted, int step) { _proposer_.update_accept_rate(is_accepted, step); }
	void reset_moves() { _proposer_.reset_moves(); }
	void adapt(int step, int num_adapt_steps) { (_proposer_.adapt)(step, num_adapt_steps); }
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

	void add_param(const std::string k, Parameter new_param);

	// Getter
	double *get_address(const std::string k) const { return params.at(k).get_address(); }
	double get_value(const std::string k) const { return params.at(k).get_value(); }
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); }

	// Setter
	void set_value(const std::string k, double new_value) { params.at(k).set_value(new_value); }
	void init_value(const std::string k, std::mt19937_64 &rng) { params.at(k).init_value(rng); }

	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); }
	void set_num_propos_pp(const std::string k) { params.at(k).set_num_propos_pp(); }
	void set_num_accept_pp(const std::string k) { params.at(k).set_num_accept_pp(); }
	void update_accept_rate(const std::string k, bool is_accepted, int step) { params.at(k).update_accept_rate(is_accepted, step); }
	void reset_moves(const std::string k) { params.at(k).reset_moves(); }

	void adapt(const std::string k, int step, int m_num_adapt_iter) { params.at(k).adapt(step, m_num_adapt_iter); }

	// Computations
	bool valid_param(const std::string k, double new_value) { return params.at(k).valid(new_value); }
	double propose(const std::string k, std::mt19937_64 &rng, double old_value) { return params.at(k).propose(rng, old_value); }
	double log_ratio_proposal(const std::string k) { return params.at(k).propose_log_ratio(); }
	double log_ratio_prior(const std::string k, double old_value, double new_value) { return params.at(k).prior_log_ratio(old_value, new_value); }

	// Display
	void print();
};

/////////////////////////////////////////////////////////
//////////// global parameters declaration //////////////
/////////////////////////////////////////////////////////

extern std::vector<double> glob_params_values;
extern Parameter_Dictionary glob_params;

/////////////////////////////////////////////////////////
//////////// Augmented : Class declaration //////////////
////////////////////////////////////////////////////////

// FIXME: to be updated

class Augmented : public Parameter
{
public:
	using Parameter::Parameter;
};

#endif // PARAMETER__H