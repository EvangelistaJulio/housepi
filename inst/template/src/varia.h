// varia.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

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

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Proposer : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME : impose mean=0, etc...

class Base_Proposer
{
protected:
	Distribution_1D *distr_ptr;
	double strat_zero;
	bool zero_start = false;
	int num_proposed_moves = 0;
	double accept_rate = 0.0;
	std::vector<int> args_to_adapt = {};

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

	double rand(sys_rnd_gen &prng) { return (distr_ptr->rand)(prng); }
	bool valid(double new_val) { return (this->*_valid_func_)(new_val); }

	virtual double propose(sys_rnd_gen &prng, double start = 0.0) { return (this->*_propose_method_)(prng, start); }
	virtual double logratio(double old_val, double new_val) { return (this->*_logratio_method_)(old_val, new_val); }

	//=====Finite Adaptive=====//
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

	std::ostringstream message_to_display();

private:
	//=====Valid Move : functional methods=====//
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
		screen::raise_error("Attempt to call a method of an undefined proposer.");
		return INF_DOUBLE;
	}

	//=====propose : functional methods=====//
	double (Base_Proposer::*_propose_method_)(sys_rnd_gen &prng, double start);

	double _undef_propose_(sys_rnd_gen &prng, double start)
	{
		UNUSED(prng);
		UNUSED(start);
		return undefined_call_error();
	}
	double _sum_to_previous_(sys_rnd_gen &prng, double start) { return start + rand(prng); }
	double _times_previous_(sys_rnd_gen &prng, double start) { return start * rand(prng); }
	double _draw_rand_(sys_rnd_gen &prng, double start)
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
	double _logratio_sum_to_zero_(double old_val, double new_val) { return _logratio_pdf_(old_val - get_strat_zero(), new_val - get_strat_zero()); }
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

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  proposer_switchboard  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Base_Proposer *proposer_switchboard(
	std::string name = UNDEFINED,
	sdargs args = {},
	double strat_zero = INF_DOUBLE,
	sdargs value_range = {});

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Undefined_Proposer : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Uniform_Proposer : class  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Uniform_Proposer : public Base_Proposer
{
public:
	//=====Constructors=====//
	Uniform_Proposer(
		sdargs new_args = {{"min", (double)-0.5}, {"max", (double)0.5}},
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {}) : Base_Proposer(new Uniform_Distr(new_args), strat_zero, value_range, {"min", "max"})
	{
		_propose_method_ = &Uniform_Proposer::_sum_to_previous_;
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
		sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {}) : Base_Proposer(new Normal_Distr(new_args), strat_zero, value_range, {1})
	{
		_propose_method_ = &Normal_Proposer::_sum_to_previous_;
		_logratio_method_ = is_zero_start() ? &Normal_Proposer::_logratio_sum_to_zero_ : &Normal_Proposer::_logratio_symm_;
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
		sdargs new_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double strat_zero = INF_DOUBLE,
		sdargs value_range = {}) : Base_Proposer(new Lognorm_Distr(new_args), strat_zero, value_range, {1})
	{
		_propose_method_ = &Lognorm_Proposer::_times_previous_;
		_logratio_method_ = is_zero_start() ? &Lognorm_Proposer::_logratio_times_zero_ : &Lognorm_Proposer::_logratio_lognorm_symm_;
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
	Prior_Proposer(Distribution_1D *distribution) : Base_Proposer(distribution)
	{
		_propose_method_ = &Prior_Proposer::_draw_rand_;
		_logratio_method_ = &Prior_Proposer::_logratio_pdf_;
	}
	~Prior_Proposer() {}
};

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Parameter : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

template <typename sdargs_int_or_ext>
class Base_Parameter
{
protected:
	double value = INF_DOUBLE;			  // value of the parameter
	bool infer = false;					  // to be inferred or static
	Distribution_1D *distr_ptr = nullptr; // prior/model distribution of the parameter
	Base_Proposer *prop_ptr = nullptr;	  // propose strategy for Monte Carlo Markov Chain

public:
	//=====Constructor=====//
	Base_Parameter(
		double value,
		std::string distr_name = UNDEFINED, sdargs_int_or_ext distr_args = {},
		std::string strat_name = UNDEFINED, sdargs strat_args = {},
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
	double distr_rand(sys_rnd_gen &prng) { return distr_ptr->rand(prng); }
	double distr_log_pdf(double new_val) { return distr_ptr->log_pdf(new_val); }
	double distr_logratio(double old_val, double new_val) { return distr_ptr->logratio(old_val, new_val); }

	//=====prop_ptr : propose methods=====//
	double propose(sys_rnd_gen &prng) { return (this->*_propose_method_)(prng); }
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
		oss << get_value() << screen::C_DELIM << ((infer) ? "true" : "false");
		if (infer)
		{
			oss << screen::C_DELIM << distr_ptr->message_to_display().str()
				<< screen::C_DELIM << prop_ptr->message_to_display().str();
		}
		return oss.str();
	}

private:
	//=====Propose : functional methods=====//
	double (Base_Parameter::*_propose_method_)(sys_rnd_gen &prng);
	double _propose_from_prior_model_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, 0); }
	double _propose_from_previous_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, get_value()); }
	double _propose_from_zero_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, get_strat_zero()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Parameter : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Parameter : public Base_Parameter<sdargs>
{
public:
	//=====Diagnostic : attributes=====//
	int dgn_sample_size = 0;	   // sample size
	long double dgn_mean = 0.;	   // mean value
	long double dgn_mean_sqr = 0.; // mean squared value

	//=====Constructor=====//
	Parameter(
		double new_val,
		std::string prior_name = UNDEFINED, sdargs prior_args = {},
		std::string strat_name = UNDEFINED, sdargs strat_args = {},
		sdargs value_range = {}) : Base_Parameter(new_val,
												  prior_name, prior_args,
												  strat_name, strat_args,
												  INF_DOUBLE, // is_zero_start() = false
												  value_range){};

	//=====Aliases=====//
	double prior_rand(sys_rnd_gen &prng) { return distr_rand(prng); }
	double prior_log_pdf(double new_val) { return distr_log_pdf(new_val); }
	double prior_logratio(double old_val, double new_val) { return distr_logratio(old_val, new_val); }

	//=====Diagnostic : methods=====//
	void init_diagnostics(int sample_size);
	void update_diagnostics();
	std::string message_to_diagnostics();
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Augmented : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

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
	double model_rand(sys_rnd_gen &prng) { return distr_rand(prng); }
	double model_log_pdf(double new_val) { return distr_log_pdf(new_val); }
	double model_logratio(double old_val, double new_val) { return distr_logratio(old_val, new_val); }
};

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Monte_Carlo_Core : class  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* This is a wrapper class containing aliases for the model params. */
class Monte_Carlo_Core
{
public:
	std::vector<std::string> param_file_required_headers = {
		"name", "value", "infer",
		"prior", "pargs",
		"strategy", "sargs", "range"};
	std::vector<std::string> param_to_be_inferred = {};
	std::vector<std::string> param_static = {};

	//=====Constructors=====//
	Monte_Carlo_Core() : _aux_unif_dist_(0.0, 1.0) {}
	virtual ~Monte_Carlo_Core(){};

	//=====Getter=====//
	void disable_assert_param();
	void assert_param(const std::string k) const;
	double get_value(const std::string k) const { return (this->*_get_value_method_)(k); }
	const double *get_address(const std::string k) const { return params.at(k).get_address(); }
	const double *get_reference() const { return &REFERENCE_VALUE; }
	bool get_infer(const std::string k) const { return params.at(k).get_infer(); }
	double get_accept_rate(const std::string k) const { return params.at(k).get_accept_rate(); }
	std::string get_message_to_display(const std::string k) { return params.at(k).message_to_display(); }
	std::string get_message_to_diagnostics(const std::string k) { return params.at(k).message_to_diagnostics(); }

	//=====Setter=====//
	void set_value(const std::string k, double new_val) { params.at(k).set_value(new_val); }
	void init_rand_params(bool is_rand_params);
	void reset_accept_rate(const std::string k) { params.at(k).reset_accept_rate(); }
	void assign_parameter(
		const std::string k, double new_val,
		std::string prior_name = UNDEFINED, sdargs prior_args = {},
		std::string strat_name = UNDEFINED, sdargs strat_args = {},
		sdargs value_range = {});
	void drop_all_params();

	//=====Parallelization=====//
	void set_num_threads(unsigned n = 1);
	int get_num_threads() const { return setup_num_threads; };

	//=====Metropolis-Hastings moves=====//
	bool propose_valid_param_move(const std::string k);
	double get_move_logratio() { return buffer_logratio; }
	void accept_param_move(const std::string k);
	void refuse_param_move(const std::string k);
	void adapt(const std::string k, double ratio) { params.at(k).adapt(ratio); }
	void init_param_diagnostics(int sample_size);
	void update_param_diagnostics();

	//=====Pseudo-Random Number Generator=====//
	void set_prng_seed(int seed);
	sys_rnd_gen &get_prng() { return (this->*_prng_method_)(); }
	double uniform() { return _aux_unif_dist_(get_prng()); }

private:
	//=====Auxiliary Attributes=====//
	std::unordered_map<std::string, Parameter> params = {}; // parameters dictionary

	double buffer_value = 0.0;								// buffer value of the current parameter
	double buffer_logratio = 0.0;							// buffer logratio (prior+proposal) of the current parameter
	std::uniform_real_distribution<double> _aux_unif_dist_; // uniform distribution
	unsigned setup_num_threads = 1;							// number of OMP threads

	//=====Getter (auxiliary methods)=====//
	double (Monte_Carlo_Core::*_get_value_method_)(const std::string k) const = &Monte_Carlo_Core::_get_value_assert_;
	double _get_value_assert_(const std::string k) const;
	double _get_value_default_(const std::string k) const;

	//=====Pseudo Random Number Generator (auxiliary methods)=====//
	sys_rnd_gen prng;							 // Mersenne Twister 19937 generator (32 or 64 bit)
	std::vector<sys_rnd_gen> prng_parallel = {}; // parallel generators (for multi-threading)
	sys_rnd_gen &(Monte_Carlo_Core::*_prng_method_)();
	sys_rnd_gen &_prng_unique_thread_() { return prng; }
#ifdef _OPENMP
	sys_rnd_gen &_prng_parallel_threads_() { return prng_parallel[get_num_threads()]; }
#endif
};

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

#endif // VARIA__H