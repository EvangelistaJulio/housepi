// varia.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#ifndef VARIA__H
#define VARIA__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>
#include <omp.h>

#include "utils.h"

// maximum number of attempts for initialization
#define MAX_NUMBER_ATTEMPTS 20

// target acceptation rate for MCMC
#define TARGET_ACCEPT_RATE 0.234

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
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {},
		std::vector<int> args_to_adapt = {});
	virtual ~Base_Proposer(){};

	//=====Getter=====//
	bool is_zero_start() const { return zero_start; }
	double get_strat_zero() const { return strat_zero; }
	double get_arg(int argID) const { return (distr_ptr->get_arg)(argID); }

	//=====Setter=====//
	void set_strat_zero(double new_val) { strat_zero = new_val; }
	void set_range(sdargs new_range);

	double rand(std::mt19937_64 &prng) { return (distr_ptr->rand)(prng); }
	bool valid(double new_val) { return (this->*_valid_func_)(new_val); }

	virtual double propose(std::mt19937_64 &prng, double start = 0.0) { return (this->*_propose_method_)(prng, start); }
	virtual double logratio(double old_val, double new_val) { return (this->*_logratio_method_)(old_val, new_val); }

	//=====Methods (for adaptive)=====//
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

	std::string message_to_display();

private:
	//=====valid : functional methods=====//
	double _range_min_ = -INF_DOUBLE;
	double _range_max_ = INF_DOUBLE;
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
		raise_error_message(err_text);
		std::exit(EXIT_FAILURE);
	}

	//=====propose : functional methods=====//
	double (Base_Proposer::*_propose_method_)(std::mt19937_64 &prng, double start);

	double _undef_propose_(std::mt19937_64 &prng, double start)
	{
		UNUSED(prng);
		UNUSED(start);
		return undefined_call_error();
	}
	double _add_to_previous_(std::mt19937_64 &prng, double start) { return start + rand(prng); }
	double _times_previous_(std::mt19937_64 &prng, double start) { return start * rand(prng); }
	double _draw_rand_(std::mt19937_64 &prng, double start)
	{
		UNUSED(start);
		return rand(prng);
	}

	//=====propose : functional methods=====//
	double (Base_Proposer::*_logratio_method_)(double old_val, double new_val);

	double _undef_logratio_(double old_val, double new_val)
	{
		UNUSED(old_val);
		UNUSED(new_val);
		return 0.0;
	}
	double _logratio_pdf_(double old_val, double new_val) { return distr_ptr->logratio(old_val, new_val); }
	// FIXME : double check these
	double _logratio_add_to_zero_(double old_val, double new_val) { return _logratio_pdf_(old_val - get_strat_zero(), new_val - get_strat_zero()); }
	double _logratio_times_zero_(double old_val, double new_val) { return _logratio_pdf_(old_val / get_strat_zero(), new_val / get_strat_zero()); }
	double _logratio_symm_(double old_val, double new_val)
	{
		UNUSED(old_val);
		UNUSED(new_val);
		return 0.0;
	}
	double _logratio_lognorm_symm_(double old_val, double new_val)
	{
		return log(old_val) - log(new_val);
	}
};

Base_Proposer *proposer_switchboard(std::string name = UNDEFINED, sdargs args = {}, double strat_zero = INF_DOUBLE, sdargs value_range = {});

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
		_logratio_method_ = &Undefined_Proposer::_undef_logratio_;
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
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {}) : Base_Proposer(new Uniform_Distr(new_args), strat_zero, value_range, {"min", "max"})
	{
		_propose_method_ = &Uniform_Proposer::_add_to_previous_;
		_logratio_method_ = &Uniform_Proposer::_logratio_symm_;
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
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {}) : Base_Proposer(new Normal_Distr(new_args), strat_zero, value_range, {1})
	{
		_propose_method_ = &Normal_Proposer::_add_to_previous_;
		_logratio_method_ = is_zero_start() ? &Normal_Proposer::_logratio_add_to_zero_ : &Normal_Proposer::_logratio_symm_;
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
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {}) : Base_Proposer(new Lognorm_Distr(new_args), strat_zero, value_range, {1})
	{
		_propose_method_ = &Lognorm_Proposer::_times_previous_;
		_logratio_method_ = is_zero_start() ? &Lognorm_Proposer::_logratio_times_zero_ : &Lognorm_Proposer::_logratio_lognorm_symm_;
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
		_logratio_method_ = &Prior_Proposer::_logratio_pdf_;
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
	double value = INF_DOUBLE;
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
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {})
	{
		// assing value
		set_value(value);

		// define distribution
		distr_ptr = distribution_switchboard(distr_name, distr_args);

		// if no strategy is defined, do not "infer"
		infer = (strat_name == UNDEFINED) ? false : true;

		// define strategy : prior/model or else with functional programming
		if (strat_name == "prior" || strat_name == "model")
		{
			prop_ptr = new Prior_Proposer(distr_ptr);
			_propose_method_ = &Base_Parameter::_propose_from_prior_model_;
		}
		else
		{
			prop_ptr = proposer_switchboard(strat_name, strat_args, strat_zero, value_range);
			_propose_method_ = is_zero_start() ? &Base_Parameter::_propose_from_zero_ : &Base_Parameter::_propose_from_previous_;
		}
	}
	virtual ~Base_Parameter(){};

	//=====Getter=====//
	double get_value() const { return value; }
	const double *get_address() const { return &value; }
	bool get_infer() const { return infer; }
	bool is_zero_start() const { return (prop_ptr->is_zero_start)(); }
	double get_strat_zero() const { return (prop_ptr->get_strat_zero)(); }

	//=====Setter=====//
	void set_value(double x) { value = x; }

	//=====distr_ptr : methods=====//
	double distr_rand(std::mt19937_64 &prng) { return distr_ptr->rand(prng); }
	double distr_log_pdf(double new_val) { return distr_ptr->log_pdf(new_val); }
	double distr_logratio(double old_val, double new_val) { return distr_ptr->logratio(old_val, new_val); }

	//=====prop_ptr : propose methods=====//
	double propose(std::mt19937_64 &prng) { return (this->*_propose_method_)(prng); }
	double propose_logratio(double old_val, double new_val) { return (prop_ptr->logratio)(old_val, new_val); }
	bool valid(double new_val) { return prop_ptr->valid(new_val); }

	//=====prop_ptr : adaptive methods=====//
	double get_accept_rate() const { return prop_ptr->get_accept_rate(); }
	void update_accept_rate(bool is_accepted) { prop_ptr->update_accept_rate(is_accepted); }
	void reset_accept_rate() { prop_ptr->reset_accept_rate(); }
	void adapt(double ratio) { (prop_ptr->adapt)(ratio); }

	//=====Display=====//
	std::string message_to_display()
	{
		std::ostringstream oss;
		oss << get_value() << DELIM << ((infer) ? "true" : "false");
		if (infer)
		{
			oss << DELIM << distr_ptr->message_to_display() << DELIM << prop_ptr->message_to_display();
		}
		return oss.str();
	}

private:
	//=====Propose : functional methods=====//
	double (Base_Parameter::*_propose_method_)(std::mt19937_64 &prng);
	double _propose_from_prior_model_(std::mt19937_64 &prng) { return (prop_ptr->propose)(prng, 0); }
	double _propose_from_previous_(std::mt19937_64 &prng) { return (prop_ptr->propose)(prng, get_value()); }
	double _propose_from_zero_(std::mt19937_64 &prng) { return (prop_ptr->propose)(prng, get_strat_zero()); }
};

/////////////////////////////////////////////////////////////////
// ///////////////////// Parameter : class /////////////////// //
/////////////////////////////////////////////////////////////////

class Parameter : public Base_Parameter<sdargs>
{
public:
	int dgn_sample_size = 0;
	long double dgn_sum = 0.;
	long double dgn_sqr_sum = 0.;

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
												  INF_DOUBLE, // is_zero_start() = false
												  value_range){};

	//=====Aliases=====//
	double prior_rand(std::mt19937_64 &prng) { return distr_rand(prng); }
	double prior_log_pdf(double new_val) { return distr_log_pdf(new_val); }
	double prior_logratio(double old_val, double new_val) { return distr_logratio(old_val, new_val); }

	void update_diagnostics();
	std::string message_to_display_diagnostics();
};

///////////////////////////////////////////////////////////
/////// Monte_Carlo_Core : class declaration //////////
///////////////////////////////////////////////////////////

// This is a wrapper class containing aliases for the model params.
class Monte_Carlo_Core
{
public:
	std::vector<std::string> param_file_required_headers = {
		"name", "value", "infer", "prior", "prior_args", "strategy", "strategy_args", "range"};
	std::vector<std::string> param_to_be_inferred = {};
	std::vector<std::string> param_static = {};

	//=====Constructors=====//
	Monte_Carlo_Core() : _aux_unif_dist_(0.0, 1.0) {}
	virtual ~Monte_Carlo_Core(){};

	// build dictionary
	void add_parameter(
		const std::string k,
		double new_val,
		std::string prior_name = UNDEFINED,
		sdargs prior_args = {},
		std::string strat_name = UNDEFINED,
		sdargs strat_args = {},
		sdargs value_range = {});

	//=====Parallelization=====//
	void set_num_threads(int n = 1);
	int get_num_threads() const { return setup_num_omp_threads; };

	//=====Pseudo-Random Number Generator=====//
	void prng_seed(int seed);
	std::mt19937_64 &get_prng() { return (this->*_prng_method_)(); }
	double uniform() { return _aux_unif_dist_(get_prng()); }

	//=====Getter=====//
	void assert_param(const std::string k);
	double get_value(const std::string k) const { return params.at(k).get_value(); }
	const double *get_address(const std::string k) const { return params.at(k).get_address(); }
	const double *get_reference() const { return &_reference; }
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); }
	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); }

	//=====Setter=====//
	void set_value(const std::string k, double new_val) { params.at(k).set_value(new_val); }
	void init_rand_params(bool is_rand_params);
	void reset_accept_rate(const std::string k) { params.at(k).reset_accept_rate(); }

	//=====Metropolis-Hastings moves=====//
	bool propose_valid_param_move(const std::string k);
	double get_move_logratio() { return buffer_logratio; }
	void accept_param_move(const std::string k);
	void refuse_param_move(const std::string k);
	void update_diagnostics();

	// acceptation rates / adaptive methods
	void adapt(const std::string k, double ratio) { params.at(k).adapt(ratio); }

	//=====Display=====//
	void display_parameters();
	void display_param_diagnostics();

private:
	std::unordered_map<std::string, Parameter> params = {};
	double buffer_value = 0.0;
	double buffer_logratio = 0.0;
	const double _reference = 1.0;

	// number of OMP threads
	int setup_num_omp_threads = 1;

	// pseudo random number generator
	std::mt19937_64 prng;
	std::vector<std::mt19937_64> prng_parallel = {};
	std::mt19937_64 &(Monte_Carlo_Core::*_prng_method_)();
	std::mt19937_64 &_prng_unique_thread_() { return prng; }
	std::mt19937_64 &_prng_parallel_threads_() { return prng_parallel[omp_get_thread_num()]; }
	std::uniform_real_distribution<double> _aux_unif_dist_;
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
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {}) : Base_Parameter(value,
												  model_name,
												  model_args_ext,
												  strat_name,
												  strat_args,
												  strat_zero,
												  value_range){};
	//=====Aliases=====//
	double model_rand(std::mt19937_64 &prng) { return distr_rand(prng); }
	double model_log_pdf(double new_val) { return distr_log_pdf(new_val); }
	double model_logratio(double old_val, double new_val) { return distr_logratio(old_val, new_val); }
};

#endif // VARIA__H