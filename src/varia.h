// varia.h //

#ifndef VARIA__H
#define VARIA__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "utils.h"

///////////////////////////////////////////////////////////////////////
// //////////////////// Base_Proposer : class ////////////////////// //
///////////////////////////////////////////////////////////////////////

// FIXME : impose mean=0, etc...

class Base_Proposer
{
protected:
	double strat_zero;
	bool zero_start = false;
	int num_proposed_moves = 0;
	double accept_rate = 0.0;
	std::vector<int> args_to_adapt = {};
	Distribution_1D *distr_ptr;

public:
	//=====Constructors=====//
	Base_Proposer(
		Distribution_1D *distr,
		double strat_zero = INFINITY,
		sdargs value_range = {},
		std::vector<int> args_to_adapt = {});
	~Base_Proposer(){};

	//=====Getter=====//
	bool is_zero_start() const { return zero_start; }
	double get_strat_zero() const { return strat_zero; }
	double get_arg(int argID) const { return (distr_ptr->get_arg)(argID); }

	//=====Setter=====//
	void set_strat_zero(double new_val) { strat_zero = new_val; }
	void set_range(sdargs new_range);

	double rand(std::mt19937_64 &rng) { return (distr_ptr->rand)(rng); }
	bool valid(double new_val) { return (this->*_valid_func_)(new_val); }

	virtual double propose(std::mt19937_64 &rng, double start = 0.0) { return (this->*_propose_method_)(rng, start); }
	virtual double log_ratio(double old_val, double new_val) { return (this->*_log_ratio_method_)(old_val, new_val); }

	//=====Methods (for adaptive)=====//
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio, double target);

	std::string display(char delim = DELIM, char dict_delim = DICT_DELIM);

private:
	//=====valid : functional methods=====//
	double _range_min_ = -INFINITY;
	double _range_max_ = INFINITY;
	bool (Base_Proposer::*_valid_func_)(double new_val);

	bool _valid_range_min_(double x) { return (x > _range_min_); }
	bool _valid_range_max_(double x) { return (x < _range_max_); }
	bool _valid_range_both_(double x) { return _valid_range_min_(x) && _valid_range_max_(x); }
	bool _valid_dummy_(double x)
	{
		UNUSED(x);
		return true;
	}

protected:
	double undefined_call_error()
	{
		std::string err_text = "Attempt to call a method of an undefined proposer.";
		exit_with_failure_error(err_text);
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
	// FIXME : double check these
	double _log_ratio_add_to_zero_(double old_val, double new_val) { return _log_ratio_pdf_(old_val - get_strat_zero(), new_val - get_strat_zero()); }
	double _log_ratio_times_zero_(double old_val, double new_val) { return _log_ratio_pdf_(old_val / get_strat_zero(), new_val / get_strat_zero()); }
	double _log_ratio_symm_(double old_val, double new_val)
	{
		UNUSED(old_val);
		UNUSED(new_val);
		return 0.0;
	}
	double _log_ratio_lognorm_symm_(double old_val, double new_val)
	{
		return log(old_val) - log(new_val);
	}
};

Base_Proposer *proposer_switchboard(std::string name = UNDEFINED, sdargs args = {}, double strat_zero = INFINITY, sdargs value_range = {});

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
		sdargs new_args = {{"min", (double)-0.5}, {"max", (double)0.5}},
		double strat_zero = INFINITY,
		sdargs value_range = {}) : Base_Proposer(new Uniform_Distr(new_args), strat_zero, value_range, {"min", "max"})
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
		sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double strat_zero = INFINITY,
		sdargs value_range = {}) : Base_Proposer(new Normal_Distr(new_args), strat_zero, value_range, {1})
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
		sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double strat_zero = INFINITY,
		sdargs value_range = {}) : Base_Proposer(new Lognorm_Distr(new_args), strat_zero, value_range, {1})
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

template <typename sdargs_int_or_ext>
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
		sdargs_int_or_ext distr_args = {},
		std::string strat_name = UNDEFINED,
		sdargs strat_args = {},
		double strat_zero = INFINITY,
		sdargs value_range = {})
	{
		set_value(value);
		distr_ptr = distribution_switchboard(distr_name, distr_args);
		infer = (strat_name == UNDEFINED) ? false : true;
		prop_ptr = (strat_name == "prior") ? (new Prior_Proposer(distr_ptr)) : proposer_switchboard(strat_name, strat_args, strat_zero, value_range);
		_propose_method_ = is_zero_start() ? &Base_Parameter::_propose_from_zero_ : &Base_Parameter::_propose_from_previous_;
	}
	~Base_Parameter() {}

	//=====Getter=====//
	double get_value() const { return value; }
	const double *get_address() const { return &value; }
	bool get_infer() const { return infer; }
	bool is_zero_start() const { return (prop_ptr->is_zero_start)(); }
	double get_strat_zero() const { return (prop_ptr->get_strat_zero)(); }

	//=====Setter=====//
	void set_value(double x) { value = x; }

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
	void adapt(double ratio, double target) { (prop_ptr->adapt)(ratio, target); }

	//=====Display=====//
	std::string display(char delim = DELIM, char dict_delim = DICT_DELIM)
	{
		std::string output = std::to_string(get_value());
		output += delim + (std::string)((infer) ? "true" : "false");
		if (infer)
		{
			output += delim + distr_ptr->display(delim, dict_delim) + delim + prop_ptr->display(delim, dict_delim);
		}
		return output;
	}

private:
	//=====Propose : functional methods=====//
	double (Base_Parameter::*_propose_method_)(std::mt19937_64 &rng);
	double _propose_from_previous_(std::mt19937_64 &rng) { return (prop_ptr->propose)(rng, get_value()); }
	double _propose_from_zero_(std::mt19937_64 &rng) { return (prop_ptr->propose)(rng, get_strat_zero()); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Parameter : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Parameter : public Base_Parameter<sdargs>
{
public:
	//=====Constructor=====//
	Parameter(
		double new_val,
		std::string prior_name = UNDEFINED,
		sdargs prior_args = {},
		std::string strat_name = UNDEFINED,
		sdargs strat_args = {},
		sdargs value_range = {}) : Base_Parameter(new_val,
												  prior_name,
												  prior_args,
												  strat_name,
												  strat_args,
												  INFINITY, // is_zero_start() = false
												  value_range){};

	//=====Aliases=====//
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
public:
	std::vector<std::string> param_file_required_headers = {
        "name", "init_value", "infer", "prior", "prior_args", "strategy", "strategy_args", "range"};
	std::vector<std::string> param_to_be_inferred = {};
	std::vector<std::string> param_static = {};
	int num_param_to_infer = 0;

	//=====Methods=====//
	void add_param(
		const std::string k,
		double new_val,
		std::string prior_name = UNDEFINED,
		sdargs prior_args = {},
		std::string strat_name = UNDEFINED,
		sdargs strat_args = {},
		sdargs value_range = {});

	//=====Getter=====//
	void check(const std::string k);
	double get_value(const std::string k) const { return params.at(k).get_value(); }
	const double *get_address(const std::string k) const { return params.at(k).get_address(); }
	const double *get_reference() const { return &_reference; }
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); }
	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); }

	//=====Setter=====//
	void set_value(const std::string k, double new_val) { params.at(k).set_value(new_val); }
	void init_rand_params(std::mt19937_64 &rng);
	void reset_accept_rate(const std::string k) { params.at(k).reset_accept_rate(); }

	//=====Metropolis-Hastings moves=====//
	bool propose_valid_move(const std::string k, std::mt19937_64 &rng);
	double get_move_log_ratio() { return buffer_log_ratio; }
	void accept_param_move(const std::string k);
	void refuse_param_move(const std::string k);

	// acceptation rates / adaptive methods
	void adapt(const std::string k, double ratio, double target) { params.at(k).adapt(ratio, target); }

	//=====Display=====//
	void display_parameters(char delim = DELIM, char dict_delim = DICT_DELIM);

private:
	std::unordered_map<std::string, Parameter> params = {};
	double buffer_value;
	double buffer_log_ratio;
	const double _reference = 1.0;
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Augmented : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Augmented : public Base_Parameter<sdargs_ext>
{
public:
	//=====Constructor=====//
	Augmented(
		double value,
		std::string model_name = UNDEFINED,
		sdargs_ext model_args_ext = {},
		std::string strat_name = UNDEFINED,
		sdargs strat_args = {},
		double strat_zero = INFINITY,
		sdargs value_range = {}) : Base_Parameter(value,
												  model_name,
												  model_args_ext,
												  strat_name,
												  strat_args,
												  strat_zero,
												  value_range){};
	//=====Aliases=====//
	double model_rand(std::mt19937_64 &rng) { return distr_rand(rng); }
	double model_log_pdf(double new_val) { return distr_log_pdf(new_val); }
	double model_log_ratio(double old_val, double new_val) { return distr_log_ratio(old_val, new_val); }
};

#endif // VARIA__H