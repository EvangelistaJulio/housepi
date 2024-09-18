/*
 * varia.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef VARIA__H
#define VARIA__H

#include <random>
#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "utils.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Proposer : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME : impose mean=0, etc...

class Base_Proposer
{
protected:
	//===== Attributes =====//
	Base_Distribution *distr_ptr;
	double from_value;

	int num_proposed_moves = 0;
	double accept_rate = 0.0;

public:
	//=====Constructors=====//
	Base_Proposer(Base_Distribution *, double rule_start_value = INF_DOUBLE,
				  sdargs range = {}, std::vector<int> args_to_adapt = {});
	virtual ~Base_Proposer() {};

	//===== Getter =====//
	std::string get_name() const { return distr_ptr->get_name(); }
	double get_arg(int argID) const { return (distr_ptr->get_arg)(argID); }
	double get_from_value() const { return from_value; }
	bool is_rule_start_from_zero() const { return _is_start_from_zero_; }

	//===== Setter =====//
	void set_rule_start_value(double x) { from_value = x; }
	void set_range(sdargs new_range);

	double rand(sys_rnd_gen &prng) { return (distr_ptr->rand)(prng); }
	bool valid(double x) const { return (this->*_valid_func_)(x); }

	virtual double propose(sys_rnd_gen &prng, double start = 0.0) { return (this->*_propose_method_)(prng, start); }
	virtual double logratio(double old_val, double new_val) const { return (this->*_logratio_method_)(old_val, new_val); }

	//=====Finite Adaptive=====//
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

	//===== Show =====//
	std::string show_arg_names() const { return distr_ptr->show_arg_names(); }
	std::string show_arg_values() const { return distr_ptr->show_arg_values(); }
	double get_range_min() const;
	double get_range_max() const;

private:
	//===== State Attributes =====//
	std::vector<int> _args_to_adapt_ = {};
	double _range_min_ = -INF_DOUBLE; // init to -infty
	double _range_max_ = INF_DOUBLE;  // init to +infty

	//=====Valid Move : functional methods=====//
	bool (Base_Proposer::*_valid_func_)(double) const = &Base_Proposer::_valid_dummy_;

	// alternative validation methods
	bool _valid_dummy_(double x) const
	{
		UNUSED(x);
		return true;
	}
	bool _valid_range_min_(double x) const { return (x > _range_min_); }
	bool _valid_range_max_(double x) const { return (x < _range_max_); }
	bool _valid_range_both_(double x) const { return _valid_range_min_(x) && _valid_range_max_(x); }

protected:
	//==== State attributes ====//
	bool _is_start_from_zero_ = false;

	//===== propose : functional methods =====//
	double (Base_Proposer::*_propose_method_)(sys_rnd_gen &, double) = &Base_Proposer::_propose_undefined_;

	// alternative proposer methods
	double _propose_undefined_(sys_rnd_gen &prng, double start)
	{
		screen::raise_error("Attempt to call a method of an undefined proposer.");
		UNUSED(prng);
		UNUSED(start);
		return std::numeric_limits<double>::quiet_NaN();
	}
	double _propose_sum_to_previous_(sys_rnd_gen &prng, double start) { return start + rand(prng); }
	double _propose_times_previous_(sys_rnd_gen &prng, double start) { return start * rand(prng); }
	double _propose_draw_rand_(sys_rnd_gen &prng, double start)
	{
		UNUSED(start);
		return rand(prng);
	}

	//===== logratio : functional methods =====//
	double (Base_Proposer::*_logratio_method_)(double, double) const = &Base_Proposer::_logratio_undefined_;

	// alternative logratios
	double _logratio_undefined_(double old_val, double new_val) const
	{
		UNUSED(old_val);
		UNUSED(new_val);
		return 0.0;
	}
	double _logratio_pdf_(double old_val, double new_val) const { return distr_ptr->logratio(old_val, new_val); }
	// FIXME : double check these
	double _logratio_sum_to_zero_(double old_val, double new_val) const { return _logratio_pdf_(old_val - get_from_value(), new_val - get_from_value()); }
	double _logratio_times_zero_(double old_val, double new_val) const { return _logratio_pdf_(old_val / get_from_value(), new_val / get_from_value()); }
	double _logratio_symm_(double old_val, double new_val) const
	{
		UNUSED(old_val);
		UNUSED(new_val);
		return 0.0;
	}
	double _logratio_lognorm_symm_(double old_val, double new_val) const { return log(old_val) - log(new_val); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  proposer_switchboard  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Base_Proposer *proposer_switchboard(
	std::string rule_name = UNDEFINED,
	sdargs rule_args = {},
	double rule_from = INF_DOUBLE,
	sdargs range = {});

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Base_Proposer
{
public:
	//=====Constructors=====//
	Uniform_Proposer(
		sdargs rule_args = {{"min", (double)-0.5}, {"max", (double)0.5}},
		double rule_from = INF_DOUBLE,
		sdargs range = {}) : Base_Proposer(new Uniform_Distribution(rule_args), rule_from, range, {"min", "max"})
	{
		_propose_method_ = &Uniform_Proposer::_propose_sum_to_previous_;
		_logratio_method_ = &Uniform_Proposer::_logratio_symm_;
	}
	~Uniform_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Normal_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Normal_Proposer : public Base_Proposer
{
public:
	//=====Constructor=====//
	Normal_Proposer(
		sdargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double rule_from = INF_DOUBLE,
		sdargs range = {}) : Base_Proposer(new Normal_Distribution(rule_args), rule_from, range, {1})
	{
		_propose_method_ = &Normal_Proposer::_propose_sum_to_previous_;
		_logratio_method_ = is_rule_start_from_zero() ? &Normal_Proposer::_logratio_sum_to_zero_ : &Normal_Proposer::_logratio_symm_;
	}
	~Normal_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Lognorm_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Lognorm_Proposer : public Base_Proposer
{
public:
	//=====Constructor=====//
	Lognorm_Proposer(
		sdargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double rule_from = INF_DOUBLE,
		sdargs range = {}) : Base_Proposer(new Lognorm_Distribution(rule_args), rule_from, range, {1})
	{
		_propose_method_ = &Lognorm_Proposer::_propose_times_previous_;
		_logratio_method_ = is_rule_start_from_zero() ? &Lognorm_Proposer::_logratio_times_zero_ : &Lognorm_Proposer::_logratio_lognorm_symm_;
	}
	~Lognorm_Proposer() {}
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Prior_Proposer : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Prior_Proposer : public Base_Proposer
{
public:
	//=====Constructor=====//
	Prior_Proposer(Base_Distribution *distribution) : Base_Proposer(distribution)
	{
		_propose_method_ = &Prior_Proposer::_propose_draw_rand_;
		_logratio_method_ = &Prior_Proposer::_logratio_pdf_;
	}
	~Prior_Proposer() {}
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Variable : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

template <typename sdargs_int_or_ptr>
class Base_Variable
{
protected:
	double value = INF_DOUBLE;				// value of the parameter
	bool infer = false;						// to be inferred or static
	Base_Distribution *distr_ptr = nullptr; // prior/model distribution of the parameter
	Base_Proposer *prop_ptr = nullptr;		// proposer for Monte Carlo Markov Chain

public:
	//=====Constructor=====//
	Base_Variable(double value, std::string distr_name = UNDEFINED, sdargs_int_or_ptr distr_args = {},
				  const std::string &rule_name = UNDEFINED, sdargs rule_args = {}, sdargs range = {}, double rule_from = INF_DOUBLE)
	{
		// assing value
		set_value(value);

		// define distribution
		distr_ptr = distribution(distr_name, distr_args);

		// if no rule is defined, do not "infer"
		infer = !screen::is_na(rule_name);

		// define proposer rule : prior/model or else with functional programming
		if (rule_name == "prior" || rule_name == "model")
		{
			prop_ptr = new Prior_Proposer(distr_ptr);
			_propose_method_ = &Base_Variable::_propose_from_prior_model_;
		}
		else
		{
			prop_ptr = proposer_switchboard(rule_name, rule_args, rule_from, range);
			_propose_method_ = is_from_zero() ? &Base_Variable::_propose_from_zero_ : &Base_Variable::_propose_from_previous_;
		}
	}
	virtual ~Base_Variable() {};

	//=====Getter=====//
	double get_value() const { return value; }
	const double *get_pointer() const { return &value; }
	bool get_infer() const { return infer; }
	bool is_from_zero() const { return (prop_ptr->is_rule_start_from_zero)(); }
	double get_from_value() const { return (prop_ptr->get_from_value)(); }

	//=====Setter=====//
	void set_value(double x) { value = x; }

	//=====distr_ptr : methods=====//
	double distr_rand(sys_rnd_gen &prng) { return distr_ptr->rand(prng); }
	double distr_log_pdf(double x) const { return distr_ptr->log_pdf(x); }
	double distr_logratio(double old_val, double new_val) const { return distr_ptr->logratio(old_val, new_val); }

	//=====prop_ptr : propose methods=====//
	double propose(sys_rnd_gen &prng) { return (this->*_propose_method_)(prng); }
	double propose_logratio(double old_val, double new_val) const { return (prop_ptr->logratio)(old_val, new_val); }
	bool valid(double x) const { return prop_ptr->valid(x); }

	//===== prop_ptr : adaptive methods =====//
	double get_accept_rate() const { return prop_ptr->get_accept_rate(); }
	void update_accept_rate(bool is_accepted) { prop_ptr->update_accept_rate(is_accepted); }
	void reset_accept_rate() { prop_ptr->reset_accept_rate(); }
	void adapt(double ratio) { (prop_ptr->adapt)(ratio); }

	//===== Show =====//
	std::string show_distr_name() const { return distr_ptr->get_name(); }
	std::string show_distr_arg_names() const { return distr_ptr->show_arg_names(); }
	std::string show_distr_arg_values() const { return distr_ptr->show_arg_values(); }
	std::string show_prop_name() const { return prop_ptr->get_name(); }
	std::string show_prop_arg_names() const { return prop_ptr->show_arg_names(); }
	std::string show_prop_arg_values() const { return prop_ptr->show_arg_values(); }
	double get_prop_range_min() const { return prop_ptr->get_range_min(); }
	double get_prop_range_max() const { return prop_ptr->get_range_max(); }

	std::string show()
	{
		using namespace screen;
		std::ostringstream oss;
		oss << get_value() << C_DELIM << truth(infer);
		if (infer)
		{
			oss << C_DELIM << show_distr_name() << C_DELIM << show_distr_arg_names() << C_DELIM << show_distr_arg_values()
				<< C_DELIM << show_prop_name() << C_DELIM << show_prop_arg_names() << C_DELIM << show_prop_arg_values()
				<< C_DELIM << get_prop_range_min() << C_DELIM << get_prop_range_max();
		}
		return oss.str();
	}

private:
	//=====Propose : functional methods=====//
	double (Base_Variable::*_propose_method_)(sys_rnd_gen &);
	double _propose_from_prior_model_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, 0); }
	double _propose_from_previous_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, get_value()); }
	double _propose_from_zero_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, get_from_value()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Parameter : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Parameter : public Base_Variable<sdargs>
{
public:
	//===== Constructor =====//
	Parameter(double new_val, const std::string &prior_name = UNDEFINED, sdargs prior_args = {},
			  const std::string &rule_name = UNDEFINED, sdargs rule_args = {},
			  sdargs range = {}) : Base_Variable(new_val, prior_name, prior_args,
												 rule_name, rule_args, range) {};

	// NOTE : default rule_from = false

	//===== Aliases =====//
	double prior_rand(sys_rnd_gen &prng) { return distr_rand(prng); }
	double prior_log_pdf(double new_val) { return distr_log_pdf(new_val); }
	double prior_logratio(double old_val, double new_val) { return distr_logratio(old_val, new_val); }

	//===== Diagnostics =====//
	void init_diagnostics(int sample_size);
	void update_diagnostics();
	std::string show_diagnostics();

private:
	//===== Diagnostics =====//
	int dgn_ss = 0;			  // sample size
	long double dgn_avg = 0.; // average value
	long double dgn_sqr = 0.; // average squared value

	long double dgn_var() { return (dgn_sqr - std::pow(dgn_avg, 2)) * dgn_ss / (dgn_ss - 1); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Augmented : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Augmented : public Base_Variable<sdargs_ptr>
{
public:
	//=====Constructor=====//
	Augmented(double value, const std::string &model_name = UNDEFINED, sdargs_ptr model_args_ptr = {},
			  const std::string &rule_name = UNDEFINED, sdargs rule_args = {}, sdargs range = {},
			  double rule_from = INF_DOUBLE) : Base_Variable(value, model_name, model_args_ptr,
															 rule_name, rule_args, range, rule_from) {};
	//=====Aliases=====//
	double model_rand(sys_rnd_gen &prng) { return distr_rand(prng); }
	double model_log_pdf(double new_val) const { return distr_log_pdf(new_val); }
	double model_logratio(double old_val, double new_val) const { return distr_logratio(old_val, new_val); }
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Monte_Carlo_Core : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* This is a wrapper class containing aliases for the model parameters. */
class Monte_Carlo_Core
{
public:
	std::vector<std::string> param_file_headers = {"name", "value", "infer",
												   "prior", "args", "argsval",
												   "rule", "rargs", "rargsval", "min", "max"};
	std::vector<std::string> param_to_be_inferred = {};
	std::vector<std::string> param_to_be_static = {};

	//===== Constructors =====//
	Monte_Carlo_Core() : _uniform_distr_(0., 1.) {}
	virtual ~Monte_Carlo_Core() {};

	//===== Getter =====//
	void disable_assert_param();
	void assert_param(const std::string k) const;
	double get_param(const std::string k) const { return (this->*_get_param_method_)(k); }
	const double *get_param_pointer(const std::string k) const { return (this->*_get_param_pointer_method_)(k); }
	const double *get_reference() const { return &REFERENCE_VALUE; }
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); }
	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); }

	//===== Setter =====//
	void set_value(const std::string k, double new_val) { params.at(k).set_value(new_val); }
	void init_rand_params(bool is_rand_params);
	void reset_accept_rate(const std::string k) { params.at(k).reset_accept_rate(); }
	void assign_parameter(const std::string k, double new_val,
						  std::string prior_name = UNDEFINED, sdargs prior_args = {},
						  std::string rule_name = UNDEFINED, sdargs rule_args = {}, sdargs range = {});
	void drop_all_params();

	//===== Parallelization =====//
	void set_num_threads(unsigned n = 1);
	int get_num_threads() const { return setup_num_threads; };

	//===== Metropolis-Hastings moves =====//
	bool propose_valid_param_move(const std::string k);
	double get_move_logratio() const { return buffer_logratio; }
	void accept_param_move(const std::string k);
	void refuse_param_move(const std::string k);
	void adapt(const std::string k, double ratio) { params.at(k).adapt(ratio); }
	void init_param_diagnostics(int sample_size);
	void update_param_diagnostics();

	//===== Pseudo-Random Number Generator =====//
	void set_prng_seed(int seed);
	sys_rnd_gen &get_prng() { return (this->*_prng_method_)(); }
	double uniform() { return _uniform_distr_(get_prng()); }

	//===== Show =====//
	std::string show_param_distr_name(const std::string k) const { return params.at(k).show_distr_name(); }
	std::string show_param_distr_arg_names(const std::string k) const { return params.at(k).show_distr_arg_names(); }
	std::string show_param_distr_arg_values(const std::string k) const { return params.at(k).show_distr_arg_values(); }
	std::string show_param_prop_name(const std::string k) const { return params.at(k).show_prop_name(); }
	std::string show_param_prop_arg_names(const std::string k) const { return params.at(k).show_prop_arg_names(); }
	std::string show_param_prop_arg_values(const std::string k) const { return params.at(k).show_prop_arg_values(); }
	double get_param_prop_range_min(const std::string k) const { return params.at(k).get_prop_range_min(); }
	double get_param_prop_range_max(const std::string k) const { return params.at(k).get_prop_range_max(); }

	std::string show_param(const std::string k) { return k + screen::C_DELIM + params.at(k).show(); }
	std::string show_param_diagnostics(const std::string k) { return k + screen::C_DELIM + params.at(k).show_diagnostics(); }

private:
	//===== Internal Attributes =====//
	std::unordered_map<std::string, Parameter> params = {}; // parameters dictionary

	double buffer_value = 0.0;								// buffer value of the current parameter
	double buffer_logratio = 0.0;							// buffer logratio (prior+proposal) of the current parameter
	std::uniform_real_distribution<double> _uniform_distr_; // uniform distribution
	unsigned setup_num_threads = 1;							// number of OMP threads

	//===== Getter (auxiliary methods) =====//
	double (Monte_Carlo_Core::*_get_param_method_)(const std::string) const = &Monte_Carlo_Core::_get_param_assert_;
	double _get_param_assert_(const std::string k) const;
	double _get_param_default_(const std::string k) const;

	const double *(Monte_Carlo_Core::*_get_param_pointer_method_)(const std::string) const = &Monte_Carlo_Core::_get_param_pointer_assert_;
	const double *_get_param_pointer_assert_(const std::string k) const;
	const double *_get_param_pointer_default_(const std::string k) const;

	//===== Pseudo Random Number Generator (auxiliary methods) =====//
	// FIXME : this should change if Rcpp according to developer guides
	sys_rnd_gen prng;							 // Mersenne Twister 19937 generator (32 or 64 bit)
	std::vector<sys_rnd_gen> prng_parallel = {}; // parallel generators (for multi-threading)
	sys_rnd_gen &(Monte_Carlo_Core::*_prng_method_)();
	sys_rnd_gen &_prng_unique_thread_() { return prng; }
#ifdef _OPENMP
	sys_rnd_gen &_prng_parallel_threads_() { return prng_parallel[omp_get_thread_num()]; }
#endif
};

#endif // VARIA__H