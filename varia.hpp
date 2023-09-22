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

#define DEFAULT_ACCEPT_RATE 0.24 // default acceptation rate for MCMC

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

	//=====Getter=====//
	bool is_zero_start() const { return zero_start; }
	double get_value_0() const { return value_0; }
	double get_arg(std::string argID) const { return (distr_ptr->get_arg)(argID); }

	//=====Setter=====//
	void set_value_0(double new_val) { value_0 = new_val; }
	void set_range(kwargs new_range);

	double rand(std::mt19937_64 &rng) { return (distr_ptr->rand)(rng); }
	bool valid(double new_val) { return (this->*_valid_func_)(new_val); }

	virtual double propose(std::mt19937_64 &rng, double start = 0.0) { return (this->*_propose_method_)(rng, start); }
	virtual double log_ratio(double old_val, double new_val) { return (this->*_log_ratio_method_)(old_val, new_val); }

	// methods for adaptive
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

private:
	//=====valid : functional methods=====//
	double _range_low_ = -INFINITY;
	double _range_up_ = INFINITY;
	bool (Base_Proposer::*_valid_func_)(double new_val);

	bool _valid_range_low_(double new_val) { return (new_val > _range_low_); }
	bool _valid_range_up_(double new_val) { return (new_val < _range_up_); }
	bool _valid_range_both_(double new_val) { return _valid_range_low_(new_val) && _valid_range_up_(new_val); }
	bool _valid_dummy_(double new_val)
	{
		UNUSED(new_val);
		return false;
	}

protected:
	double undefined_call_error()
	{
		std::cerr << "\nERROR: Attempt to call a method of an undefined proposer." << std::endl;
		std::cerr << "Exit...\n"
				  << std::endl;
		std::exit(EXIT_FAILURE);
		return INFINITY;
	}

	//=====propose : functional methods=====//
	double (Base_Proposer::*_propose_method_)(std::mt19937_64 &rng, double start);

	double _undef_propose_(std::mt19937_64 &rng, double start)
	{
		UNUSED(rng);
		UNUSED(start);
		return undefined_call_error();
	}
	double _add_to_previous_(std::mt19937_64 &rng, double start) { return start + rand(rng); }
	double _times_previous_(std::mt19937_64 &rng, double start) { return start * rand(rng); }
	double _draw_rand_(std::mt19937_64 &rng, double start)
	{
		UNUSED(start);
		return rand(rng);
	}

	//=====propose : functional methods=====//
	double (Base_Proposer::*_log_ratio_method_)(double old_val, double new_val);

	double _undef_log_ratio_(double old_val, double new_val)
	{
		UNUSED(old_val);
		UNUSED(new_val);
		return undefined_call_error();
	}
	double _log_ratio_pdf_(double old_val, double new_val) { return distr_ptr->log_ratio(old_val, new_val); }
	double _log_ratio_add_to_zero_(double old_val, double new_val) { return _log_ratio_pdf_(old_val - get_value_0(), new_val - get_value_0()); }
	double _log_ratio_times_zero_(double old_val, double new_val) { return _log_ratio_pdf_(old_val / get_value_0(), new_val / get_value_0()); }
	double _log_ratio_symm_(double old_val, double new_val)
	{
		UNUSED(old_val);
		UNUSED(new_val);
		return 0.0;
	}
	double _log_ratio_lognorm_symm_(double old_val, double new_val)
	{
		return (std::log(old_val) - std::log(new_val)); // FIXME : check *2
	}
};

Base_Proposer *proposer_switchboard(std::string name = UNDEFINED, kwargs args = {}, double value_0 = INFINITY, kwargs value_range = {});

///////////////////////////////////////////////////////////////////////
// /////////////////// Undefined_Proposer : class ////////////////// //
///////////////////////////////////////////////////////////////////////

class Undefined_Proposer : public Base_Proposer
{
public:
	//=====Constructors=====//
	Undefined_Proposer() : Base_Proposer(new Undefined_Distr())
	{
		_propose_method_ = &Undefined_Proposer::_undef_propose_;
		_log_ratio_method_ = &Undefined_Proposer::_undef_log_ratio_;
	}
	~Undefined_Proposer() {}
};

///////////////////////////////////////////////////////////////////////
// ///////////////////// Uniform_Proposer : class /////////////////// //
///////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Base_Proposer
{
public:
	//=====Constructors=====//
	Uniform_Proposer(
		kwargs new_args = {{"low", (double)-0.5}, {"up", (double)0.5}},
		double value_0 = INFINITY,
		kwargs value_range = {}) : Base_Proposer(new Uniform_Distr(new_args), value_0, value_range, {"low", "up"})
	{
		_propose_method_ = &Uniform_Proposer::_add_to_previous_;
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
	//=====Constructor=====//
	Normal_Proposer(
		kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double value_0 = INFINITY,
		kwargs value_range = {}) : Base_Proposer(new Normal_Distr(new_args), value_0, value_range, {"sd"})
	{
		_propose_method_ = &Normal_Proposer::_add_to_previous_;
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
	//=====Constructor=====//
	Lognorm_Proposer(
		kwargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double value_0 = INFINITY,
		kwargs value_range = {}) : Base_Proposer(new Lognorm_Distr(new_args), value_0, value_range, {"sd"})
	{
		_propose_method_ = &Lognorm_Proposer::_times_previous_;
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
	//=====Constructor=====//
	Prior_Proposer(Distribution_1D *distribution) : Base_Proposer(distribution)
	{
		_propose_method_ = &Prior_Proposer::_draw_rand_;
		_log_ratio_method_ = &Prior_Proposer::_log_ratio_pdf_;
	}
	~Prior_Proposer() {}
};

///////////////////////////////////////////////////////////////////////
// ///////////////// Base_Parameter  : class /////////////////// //
///////////////////////////////////////////////////////////////////////

template <typename kwargs_int_or_ext>
class Base_Parameter
{
protected:
	double value = INFINITY;
	bool infer = false;
	Distribution_1D *distr_ptr = nullptr;
	Base_Proposer *prop_ptr = nullptr;

public:
	//=====Constructor=====//
	Base_Parameter(
		double value,
		std::string distr_name = UNDEFINED,
		kwargs_int_or_ext distr_args = {},
		std::string strat_name = UNDEFINED,
		kwargs strat_args = {},
		double value_0 = INFINITY,
		kwargs value_range = {})
	{
		set_value(value);
		distr_ptr = distribution_switchboard(distr_name, distr_args);
		infer = (strat_name == UNDEFINED) ? false : true;
		prop_ptr = (strat_name == "prior") ? (new Prior_Proposer(distr_ptr)) : proposer_switchboard(strat_name, strat_args, value_0, value_range);
		_propose_method_ = is_zero_start() ? &Base_Parameter::_propose_from_zero_ : &Base_Parameter::_propose_from_previous_;
	}
	~Base_Parameter() {}

	//=====Getter=====//
	double get_value() const { return value; }
	const double *get_address() const { return &value; }
	bool get_infer() const { return infer; }
	bool is_zero_start() const { return (prop_ptr->is_zero_start)(); }
	double get_value_0() const { return (prop_ptr->get_value_0)(); }

	//=====Setter=====//
	void set_value(double value) { this->value = value; }

	//=====distr_ptr : methods=====//
	double distr_rand(std::mt19937_64 &rng) { return distr_ptr->rand(rng); }
	double distr_log_pdf(double new_val) { return distr_ptr->log_pdf(new_val); }
	double distr_log_ratio(double old_val, double new_val) { return distr_ptr->log_ratio(old_val, new_val); }

	//=====prop_ptr : propose methods=====//
	double propose(std::mt19937_64 &rng) { return (this->*_propose_method_)(rng); }
	double propose_log_ratio(double old_val, double new_val) { return (prop_ptr->log_ratio)(old_val, new_val); }
	bool valid(double new_val) { return prop_ptr->valid(new_val); }

	//=====prop_ptr : adaptive methods=====//
	double get_accept_rate() const { return prop_ptr->get_accept_rate(); }
	void update_accept_rate(bool is_accepted) { prop_ptr->update_accept_rate(is_accepted); }
	void reset_accept_rate() { prop_ptr->reset_accept_rate(); }
	void adapt(double ratio) { (prop_ptr->adapt)(ratio); }

private:
	//=====Propose : functional methods=====//
	double (Base_Parameter::*_propose_method_)(std::mt19937_64 &rng);
	double _propose_from_previous_(std::mt19937_64 &rng) { return (prop_ptr->propose)(rng, get_value()); }
	double _propose_from_zero_(std::mt19937_64 &rng) { return (prop_ptr->propose)(rng, get_value_0()); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Parameter : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Parameter : public Base_Parameter<kwargs>
{
public:
	//=====Constructor=====//
	Parameter(
		double new_val,
		std::string prior_name = UNDEFINED,
		kwargs prior_args = {},
		std::string strat_name = UNDEFINED,
		kwargs strat_args = {},
		kwargs value_range = {}) : Base_Parameter(new_val,
												  prior_name,
												  prior_args,
												  strat_name,
												  strat_args,
												  INFINITY, // is_zero_start() = false
												  value_range){};

	// Aliases
	double prior_rand(std::mt19937_64 &rng) { return distr_rand(rng); }
	double prior_log_pdf(double new_val) { return distr_log_pdf(new_val); }
	double prior_log_ratio(double old_val, double new_val) { return distr_log_ratio(old_val, new_val); }
};

///////////////////////////////////////////////////////////
/////// Parameter_Dictionary : class declaration //////////
///////////////////////////////////////////////////////////

// This is a wrapper class containing aliases for the model params.
class Parameter_Dictionary
{
private:
	std::unordered_map<std::string, Parameter> params = {};
	double buffer_value;
	double buffer_log_ratio;

public:
	std::vector<std::string> param_to_be_inferred = {};
	std::vector<std::string> param_not_to_be_inferred = {};

	//=====Methods=====//
	void add_param(
		const std::string k,
		double new_val,
		std::string prior_name = UNDEFINED,
		kwargs prior_args = {},
		std::string strat_name = UNDEFINED,
		kwargs strat_args = {},
		kwargs value_range = {});

	//=====Getter=====//
	double get_value(const std::string k) const { return params.at(k).get_value(); }
	const double *get_address(const std::string k) const { return params.at(k).get_address(); }
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); }
	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); }

	//=====Setter=====//
	void set_value(const std::string k, double new_val) { params.at(k).set_value(new_val); }
	void init_rand_params(std::mt19937_64 &rng);
	void reset_accept_rate(const std::string k) { params.at(k).reset_accept_rate(); }

	//=====Metropolis-Hastings moves=====//
	bool propose_valid_move(const std::string k, std::mt19937_64 &rng);
	double get_move_log_ratio() { return buffer_log_ratio; }
	void accept_move(const std::string k);
	void refuse_move(const std::string k);

	// acceptation rates / adaptive methods
	void adapt(const std::string k, double ratio) { params.at(k).adapt(ratio); }

	//=====Display=====//
	void display_parameters(char delim = ' ');
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Augmented : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Augmented : public Base_Parameter<kwargs_ext>
{
public:
	//=====Constructor=====//
	Augmented(
		double value,
		std::string model_name = UNDEFINED,
		kwargs_ext model_args_ext = {},
		std::string strat_name = UNDEFINED,
		kwargs strat_args = {},
		double value_0 = INFINITY,
		kwargs value_range = {}) : Base_Parameter(value,
												  model_name,
												  model_args_ext,
												  strat_name,
												  strat_args,
												  value_0,
												  value_range){};
	//=====Aliases=====//
	double model_rand(std::mt19937_64 &rng) { return distr_rand(rng); }
	double model_log_pdf(double new_val) { return distr_log_pdf(new_val); }
	double model_log_ratio(double old_val, double new_val) { return distr_log_ratio(old_val, new_val); }
};

#endif // VARIA__H