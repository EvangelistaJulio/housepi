// parameter.hpp //

#ifndef PARAMETER__H
#define PARAMETER__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "utils.hpp"

///////////////////////////////////////////////////////////////////////
// //////////////////// Base_Proposer : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

// FIXME : certain argument must be equal to 0.

class Base_Proposer
{
protected:
	// acceptation
	int num_accept = 0;
	int num_propos = 0;
	double accept_rate = 0.0;

	double _adapter_(int step, int num_adapt_iter);

public:
	// Constructors
	Base_Proposer() = default;

	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted, int step);
	void set_num_accept_pp() { num_accept++; };
	void set_num_propos_pp() { num_propos++; };
	void reset_moves();

	void set_range(kwargs new_range);
	bool valid(double new_value) { return (this->*_valid_func_)(new_value); }

	// --- functions (virtual) --- //
	virtual double propose(std::mt19937_64 &rng, double old_value) = 0;
	virtual double log_ratio() = 0; // NOTE : the options at the moment do not depend on old/new values
	virtual void adapt(int step, int num_adapt_iter) = 0;

private:
	// --- Validation : functional methods --- //
	double _range_low_ = -INFINITY;
	double _range_up_ = INFINITY;
	bool (Base_Proposer::*_valid_func_)(double new_value);

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
// ///////////////////// Normal_Proposer : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Base_Proposer
{
private:
	Normal_Distr distrib;
	std::vector<std::string> args_to_adapt = {"sd"};

public:
	// Constructors
	Normal_Proposer(
		kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		kwargs value_range = {}) : distrib(new_args) { set_range(value_range); }
	~Normal_Proposer() {}

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
// ///////////////////// Lognorm_Proposer : class ////////////////// //
///////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Base_Proposer
{
private:
	Lognorm_Distr distrib;
	std::vector<std::string> args_to_adapt = {"sd"};

public:
	// Constructors
	Lognorm_Proposer(
		kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		kwargs value_range = {}) : distrib(new_args) { set_range(value_range); }
	~Lognorm_Proposer() {}

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
// /////////////////// Base_Parameter : class ///////////////// //
/////////////////////////////////////////////////////////////////

class Base_Parameter
{
private:
	// objects for inference
	Distribution_1D *distribution_instance;
	Base_Proposer *proposer_instance;

public:
	// Constructor // UPDATE: templates are messy in header + constructor to be cleaned
	Base_Parameter(
		Distribution_1D *distribution,
		Base_Proposer *proposer) : distribution_instance(distribution), proposer_instance(proposer) {}

	// Computations
	bool valid(double new_value) { return proposer_instance->valid(new_value); }

	// Methods for distribution_instance
	double distr_rand(std::mt19937_64 &rng) { return distribution_instance->rand(rng); }
	double distr_log_pdf(double new_value) { return distribution_instance->log_pdf(new_value); }
	double distr_log_ratio(double old_value, double new_value) { return distribution_instance->log_ratio(old_value, new_value); }

	// Methods for proposer_instance
	double propose(std::mt19937_64 &rng, double old_value) { return (proposer_instance->propose)(rng, old_value); }
	double propose_log_ratio() { return (proposer_instance->log_ratio)(); }

	double get_accept_rate() const { return proposer_instance->get_accept_rate(); }
	void set_num_propos_pp() { proposer_instance->set_num_propos_pp(); }
	void set_num_accept_pp() { proposer_instance->set_num_accept_pp(); }
	void update_accept_rate(bool is_accepted, int step) { proposer_instance->update_accept_rate(is_accepted, step); }
	void reset_moves() { proposer_instance->reset_moves(); }
	void adapt(int step, int num_adapt_steps) { (proposer_instance->adapt)(step, num_adapt_steps); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Parameter : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Parameter
{
public:
	Base_Parameter _base_;
	double *p_value_ptr;
	bool p_infer = false;

	Parameter(
		double *new_value_ptr,
		std::string prior_name = UNDEF,
		kwargs prior_args = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		kwargs value_range = {});

	// Getter
	bool get_infer() const { return p_infer; }
	double get_value() const { return *(p_value_ptr); }
	double *get_address() const { return p_value_ptr; }

	// Setter
	void set_value(double new_value) { *(p_value_ptr) = new_value; }

	// Aliases
	double prior_rand(std::mt19937_64 &rng) { return _base_.distr_rand(rng); }
	double prior_log_ratio(double old_value, double new_value) { return _base_.distr_log_ratio(old_value, new_value); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Augmented : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Augmented
{
public:
	Base_Parameter _base_;
	double a_value;
	double a_value_0;
	bool a_infer = false;

	Augmented(
		double new_value,
		std::string model_name = UNDEF,
		kwargs_ext model_args_ext = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		kwargs value_range = {});

	// Getter
	bool get_infer() const { return a_infer; }
	double get_value() const { return a_value; }
	double get_value_0() const { return a_value_0; }

	// Setter
	void set_value(double new_value) { a_value = new_value; }

	// Aliases
	double model_rand(std::mt19937_64 &rng) { return _base_.distr_rand(rng); }
	double model_log_pdf(double new_value) { return _base_.distr_log_pdf(new_value); }
	double model_log_ratio(double old_value, double new_value) { return _base_.distr_log_ratio(old_value, new_value); }
};

///////////////////////////////////////////////////////////
/////// Parameter_Dictionary : class declaration //////////
///////////////////////////////////////////////////////////

// This is a wrapper class containing aliases for the global params.
class Parameter_Dictionary
{
private:
	// Attributes
	std::vector<double> param_values = {};
	std::unordered_map<std::string, Parameter> params = {};
	size_t _counter_ = 0;

public:
	// Attributes
	std::vector<std::string> to_be_inferred = {};
	std::vector<std::string> not_to_be_inferred = {};

	void resize(size_t num) { param_values.resize(num); }
	void add_param(
		const std::string k,
		double new_value,
		std::string prior_name = UNDEF,
		kwargs prior_args = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		kwargs value_range = {});

	// Getter
	double *get_address(const std::string k) const { return params.at(k).get_address(); }
	double get_value(const std::string k) const { return params.at(k).get_value(); }
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); }

	// Setter
	void set_value(const std::string k, double new_value) { params.at(k).set_value(new_value); }
	void init_value(const std::string k, std::mt19937_64 &rng) { set_value(k, params.at(k).prior_rand(rng)); }

	double get_accept_rate(const std::string k) const { return params.at(k)._base_.get_accept_rate(); }
	void set_num_propos_pp(const std::string k) { params.at(k)._base_.set_num_propos_pp(); }
	void set_num_accept_pp(const std::string k) { params.at(k)._base_.set_num_accept_pp(); }
	void update_accept_rate(const std::string k, bool is_accepted, int step) { params.at(k)._base_.update_accept_rate(is_accepted, step); }
	void reset_moves(const std::string k) { params.at(k)._base_.reset_moves(); }

	void adapt(const std::string k, int step, int m_num_adapt_iter) { params.at(k)._base_.adapt(step, m_num_adapt_iter); }

	// Computations
	bool valid_param(const std::string k, double new_value) { return params.at(k)._base_.valid(new_value); }
	double propose(const std::string k, std::mt19937_64 &rng, double old_value) { return params.at(k)._base_.propose(rng, old_value); }
	double log_ratio_proposal(const std::string k) { return params.at(k)._base_.propose_log_ratio(); }
	double log_ratio_prior(const std::string k, double old_value, double new_value) { return params.at(k).prior_log_ratio(old_value, new_value); }

	// Display
	void print();
};

/////////////////////////////////////////////////////////
//////////// global parameters declaration //////////////
/////////////////////////////////////////////////////////

extern Parameter_Dictionary glob_params;

#endif // PARAMETER__H