// varia.hpp //

#ifndef VARIA__H
#define VARIA__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "utils.hpp"

#define DEFAULT_ACCEPT_RATE 0.25 // default acceptation rate for MCMC

///////////////////////////////////////////////////////////////////////
// //////////////////// Base_Proposer : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

// FIXME : impose mean=0, etc...

class Base_Proposer
{
protected:
	double value_0;
	bool zero_start = false;
	unsigned num_proposed_moves = 0;
	double accept_rate = 0.0;
	std::vector<std::string> args_to_adapt = {};
	Distribution_1D *distr_ptr;

public:
	// Constructors
	Base_Proposer(
		Distribution_1D *distr,
		double value_0 = INFINITY,
		kwargs value_range = {},
		std::vector<std::string> args_to_adapt = {});
	~Base_Proposer(){};

	// --- Getter --- //
	bool is_zero_start() const { return zero_start; }
	double get_value_0() const { return value_0; }
	double get_arg(std::string argID) const { return (distr_ptr->get_arg)(argID); }

	// --- Setter --- //
	void set_value_0(double new_value) { value_0 = new_value; }
	void set_range(kwargs new_range);

	double rand(std::mt19937_64 &rng) { return (distr_ptr->rand)(rng); }
	bool invalid(double new_value) { return (this->*_invalid_func_)(new_value); }

	virtual double propose(std::mt19937_64 &rng, double start_value = 0.0) { return (this->*_propose_method_)(rng, start_value); }
	virtual double log_ratio(double old_value, double new_value) { return (this->*_log_ratio_method_)(old_value, new_value); }

	// methods for adaptive
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

private:
	// --- invalid : functional methods --- //
	double _range_low_ = -INFINITY;
	double _range_up_ = INFINITY;
	bool (Base_Proposer::*_invalid_func_)(double new_value);

	bool _invalid_range_low_(double new_value) { return (new_value <= _range_low_); }
	bool _invalid_range_up_(double new_value) { return (new_value >= _range_up_); }
	bool _invalid_range_both_(double new_value) { return _invalid_range_low_(new_value) && _invalid_range_up_(new_value); }
	bool _invalid_dummy_(double new_value)
	{
		(void)new_value;
		return false;
	}

protected:
	// --- propose : functional methods --- //
	double (Base_Proposer::*_propose_method_)(std::mt19937_64 &rng, double start_value);
	double _draw_rand_(std::mt19937_64 &rng, double start_value)
	{
		(void)start_value;
		return rand(rng);
	}
	double _add_to_zero_(std::mt19937_64 &rng, double start_value)
	{
		(void)start_value;
		return get_value_0() + rand(rng);
	}
	double _add_to_previous_(std::mt19937_64 &rng, double start_value)
	{
		return start_value + rand(rng);
	}
	double _times_zero_(std::mt19937_64 &rng, double start_value)
	{
		(void)start_value;
		return get_value_0() * rand(rng);
	}
	double _times_previous_(std::mt19937_64 &rng, double start_value)
	{
		return start_value * rand(rng);
	}

	// --- propose : functional methods --- //
	double (Base_Proposer::*_log_ratio_method_)(double old_value, double new_value);

	double _log_ratio_pdf_(double old_value, double new_value)
	{
		return distr_ptr->log_ratio(old_value, new_value);
	}
	double _log_ratio_add_to_zero_(double old_value, double new_value)
	{
		return _log_ratio_pdf_(old_value - get_value_0(), new_value - get_value_0());
	}
	double _log_ratio_times_zero_(double old_value, double new_value)
	{
		return _log_ratio_pdf_(old_value / get_value_0(), new_value / get_value_0());
	}
	double _log_ratio_symm_(double old_value, double new_value)
	{
		(void)old_value;
		(void)new_value;
		return 0.0;
	}
	double _log_ratio_lognorm_symm_(double old_value, double new_value)
	{
		return (std::log(old_value) - std::log(new_value)); // FIXME : check *2
	}
};

Base_Proposer *proposer_switchboard(std::string name, kwargs args = {}, double value_0 = INFINITY, kwargs value_range = {});

///////////////////////////////////////////////////////////////////////
// ///////////////////// Uniform_Proposer : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Base_Proposer
{
public:
	// Constructors
	Uniform_Proposer(
		kwargs new_args = {{"low", (double)-0.5}, {"up", (double)0.5}},
		double value_0 = INFINITY,
		kwargs value_range = {}) : Base_Proposer(new Uniform_Distr(new_args), value_0, value_range, {"low", "up"})
	{
		_propose_method_ = is_zero_start() ? &Uniform_Proposer::_add_to_zero_ : &Uniform_Proposer::_add_to_previous_;
		_log_ratio_method_ = &Uniform_Proposer::_log_ratio_symm_;
	}
	~Uniform_Proposer() {}
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Normal_Proposer : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Base_Proposer
{
public:
	// Constructors
	Normal_Proposer(
		kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double value_0 = INFINITY,
		kwargs value_range = {}) : Base_Proposer(new Normal_Distr(new_args), value_0, value_range, {"sd"})
	{
		_propose_method_ = is_zero_start() ? &Normal_Proposer::_add_to_zero_ : &Normal_Proposer::_add_to_previous_;
		_log_ratio_method_ = is_zero_start() ? &Normal_Proposer::_log_ratio_add_to_zero_ : &Normal_Proposer::_log_ratio_symm_;
	}
	~Normal_Proposer() {}
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Lognorm_Proposer : class ////////////////// //
///////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Base_Proposer
{
public:
	// Constructors
	Lognorm_Proposer(
		kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double value_0 = INFINITY,
		kwargs value_range = {}) : Base_Proposer(new Lognorm_Distr(new_args), value_0, value_range, {"sd"})
	{
		_propose_method_ = is_zero_start() ? &Lognorm_Proposer::_times_zero_ : &Lognorm_Proposer::_times_previous_;
		_log_ratio_method_ = is_zero_start() ? &Lognorm_Proposer::_log_ratio_times_zero_ : &Lognorm_Proposer::_log_ratio_lognorm_symm_;
	}
	~Lognorm_Proposer() {}
};

///////////////////////////////////////////////////////////////////////
// ////////////////////// Prior_Proposer : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Prior_Proposer : public Base_Proposer
{
public:
	// Constructors
	Prior_Proposer(Distribution_1D *distribution) : Base_Proposer(distribution)
	{
		_propose_method_ = &Prior_Proposer::_draw_rand_;
		_log_ratio_method_ = &Prior_Proposer::_log_ratio_pdf_;
	}
	~Prior_Proposer() {}
};

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////
// ////////////////// Base_Parameter : class ///////////////// //
/////////////////////////////////////////////////////////////////

class Base_Parameter
{
public:
	// objects for inference
	Distribution_1D *distr_ptr;
	Base_Proposer *prop_ptr;

	// Constructors
	Base_Parameter(Distribution_1D *distr) : distr_ptr(distr), prop_ptr(new Prior_Proposer(distr_ptr)) {}
	Base_Parameter(Distribution_1D *distr, Base_Proposer *prop) : distr_ptr(distr), prop_ptr(prop) {}
};

/////////////////////////////////////////////////////////////////
// ////////////// Wrap_Base_Parameter : class //////////////// //
/////////////////////////////////////////////////////////////////

template <typename kw_doub_ptr>
class Wrap_Base_Parameter
{
protected:
	Base_Parameter _base_;
	bool infer = false;

public:
	Wrap_Base_Parameter(
		std::string prior_name = UNDEF,
		kw_doub_ptr prior_args = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		double value_0 = INFINITY,
		kwargs value_range = {}) : _base_(nullptr, nullptr)
	{
		// --- metropolis strategy --- //

		if (strat_name != UNDEF)
		{
			Base_Proposer *proposer;
			Distribution_1D *distribution;

			infer = true;
			distribution = distribution_switchboard(prior_name, prior_args);

			// --- proposer ---//
			if (strat_name == "prior")
			{
				proposer = new Prior_Proposer(distribution);
			}
			else
			{
				proposer = proposer_switchboard(strat_name, strat_args, value_0, value_range);
			}
			_base_ = Base_Parameter(distribution, proposer);
		}
	}

	// --- Methods --- //

	// --- Getter --- //
	bool get_infer() const { return infer; }
	bool is_zero_start() const { return (_base_.prop_ptr->is_zero_start)(); }
	double get_value_0() const { return (_base_.prop_ptr->get_value_0)();}

	// --- distr_ptr : methods --- //
	double distr_rand(std::mt19937_64 &rng) { return _base_.distr_ptr->rand(rng); }
	double distr_log_pdf(double new_value) { return _base_.distr_ptr->log_pdf(new_value); }
	double distr_log_ratio(double old_value, double new_value) { return _base_.distr_ptr->log_ratio(old_value, new_value); }

	// --- prop_ptr : methods --- //
	double propose_from(std::mt19937_64 &rng, double start_value) { return (_base_.prop_ptr->propose)(rng, start_value); }
	double propose_log_ratio(double old_value, double new_value) { return (_base_.prop_ptr->log_ratio)(old_value, new_value); }
	bool invalid(double new_value) { return _base_.prop_ptr->invalid(new_value); }
	double get_accept_rate() const { return _base_.prop_ptr->get_accept_rate(); }
	void update_accept_rate(bool is_accepted) { _base_.prop_ptr->update_accept_rate(is_accepted); }
	void reset_accept_rate() { _base_.prop_ptr->reset_accept_rate(); }
	void adapt(double ratio) { (_base_.prop_ptr->adapt)(ratio); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Parameter : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Parameter : public Wrap_Base_Parameter<kwargs>
{
private:
	double *p_value_ptr;

public:
	// --- Constructor --- //
	Parameter(
		double *new_value_ptr,
		std::string prior_name = UNDEF,
		kwargs prior_args = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		kwargs value_range = {}) : Wrap_Base_Parameter(prior_name,
													   prior_args,
													   strat_name,
													   strat_args,
													   INFINITY, // is_zero_start() = false
													   value_range)
	{
		p_value_ptr = new_value_ptr;
		if (get_infer() == true)
		{
			_propose_method_ = is_zero_start() ? &Parameter::_propose_from_zero_ : &Parameter::_propose_from_previous_;
		}
	};

	// --- Getter --- //
	double get_value() const { return *(p_value_ptr); }
	double *get_address() const { return p_value_ptr; }

	// Setter
	void set_value(double new_value) { *(p_value_ptr) = new_value; }

	// Aliases
	double prior_rand(std::mt19937_64 &rng) { return distr_rand(rng); }
	double prior_log_pdf(double new_value) { return distr_log_pdf(new_value); }
	double prior_log_ratio(double old_value, double new_value) { return distr_log_ratio(old_value, new_value); }
	double propose(std::mt19937_64 &rng) { return (this->*_propose_method_)(rng); }

private:
	// --- Propose : functional methods --- //
	double (Parameter::*_propose_method_)(std::mt19937_64 &rng);
	double _propose_from_previous_(std::mt19937_64 &rng) { return propose_from(rng, get_value()); }
	double _propose_from_zero_(std::mt19937_64 &rng) { return propose_from(rng, get_value_0()); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Augmented : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Augmented : public Wrap_Base_Parameter<kwargs_ext>
{
private:
	double a_value;

public:
	// --- Constructor --- //
	Augmented(
		double new_value,
		std::string model_name = UNDEF,
		kwargs_ext model_args_ext = {},
		std::string strat_name = UNDEF,
		kwargs strat_args = {},
		double value_0 = INFINITY,
		kwargs value_range = {}) : Wrap_Base_Parameter(model_name,
													   model_args_ext,
													   strat_name,
													   strat_args,
													   value_0,
													   value_range)
	{

		a_value = new_value;
		if (get_infer() == true)
		{
			_propose_method_ = is_zero_start() ? &Augmented::_propose_from_zero_ : &Augmented::_propose_from_previous_;
		}
	};

	// --- Getter --- //
	double get_value() const { return a_value; }

	// Setter
	void set_value(double new_value) { a_value = new_value; }

	// Aliases
	double model_rand(std::mt19937_64 &rng) { return distr_rand(rng); }
	double model_log_pdf(double new_value) { return distr_log_pdf(new_value); }
	double model_log_ratio(double old_value, double new_value) { return distr_log_ratio(old_value, new_value); }
	double propose(std::mt19937_64 &rng) { return (this->*_propose_method_)(rng); }

private:
	// --- Propose : functional methods --- //
	double (Augmented::*_propose_method_)(std::mt19937_64 &rng);
	double _propose_from_previous_(std::mt19937_64 &rng) { return propose_from(rng, get_value()); }
	double _propose_from_zero_(std::mt19937_64 &rng) { return propose_from(rng, get_value_0()); }
};

///////////////////////////////////////////////////////////
/////// Parameter_Dictionary : class declaration //////////
///////////////////////////////////////////////////////////

// This is a wrapper class containing aliases for the global params.
class Parameter_Dictionary
{
private:
	std::vector<double> param_values = {}; // FIXME can I move this inside Parameter?
	std::unordered_map<std::string, Parameter> params = {};
	size_t _counter_ = 0;

public:
	std::vector<std::string> to_be_inferred = {};
	std::vector<std::string> not_to_be_inferred = {};

	// --- Methods --- //

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
	void random_init_values(std::mt19937_64 &rng);

	// acceptation rates / adaptive methods
	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); }
	void update_accept_rate(const std::string k, bool is_accepted) { params.at(k).update_accept_rate(is_accepted); }
	void reset_accept_rate(const std::string k) { params.at(k).reset_accept_rate(); }
	void adapt(const std::string k, double ratio) { params.at(k).adapt(ratio); }

	// Computations
	bool invalid_param(const std::string k, double new_value) { return params.at(k).invalid(new_value); }
	double propose(const std::string k, std::mt19937_64 &rng) { return params.at(k).propose(rng); }
	double log_ratio_proposal(const std::string k, double old_value, double new_value) { return params.at(k).propose_log_ratio(old_value, new_value); }
	double log_ratio_prior(const std::string k, double old_value, double new_value) { return params.at(k).prior_log_ratio(old_value, new_value); }

	// Display
	void print(char delim = ' ');
};

/////////////////////////////////////////////////////////
//////////// global parameters declaration //////////////
/////////////////////////////////////////////////////////

extern Parameter_Dictionary glob_params;

#endif // VARIA__H