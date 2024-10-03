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
#include <map>
#include <unordered_map>
#include <cmath>

#include "distributions.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

// bool Invalid_Proposer_Range_Warning = false; // warning for invalid proposer range declaration

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
public:
	//=====Constructors=====//
	Base_Proposer(Base_Distribution *, double rule_start_value = INF_DOUBLE,
				  vector<int> args_to_adapt = {});
	virtual ~Base_Proposer() {};

	//===== Getter =====//
	string get_name() const { return distr_ptr->get_name(); }
	double get_arg(int argID) const { return (distr_ptr->get_arg)(argID); }
	double get_from_value() const { return from_value; }
	bool is_rule_start_from_zero() const { return _is_start_from_zero_; }

	//===== Setter =====//
	void set_rule_start_value(double x) { from_value = x; }

	double rand(sys_rnd_gen &prng) { return (distr_ptr->rand)(prng); }

	virtual double propose(sys_rnd_gen &prng, double start = 0.0) { return (this->*_propose_method_)(prng, start); }
	virtual double logratio(double old_val, double new_val) const { return (this->*_logratio_method_)(old_val, new_val); }

	//=====Finite Adaptive=====//
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

	//===== Show =====//
	string show_arg_names() const { return distr_ptr->show_arg_names(); }
	string show_arg_values() const { return distr_ptr->show_arg_values(); }

protected:
	//===== Attributes =====//
	Base_Distribution *distr_ptr;
	double from_value = INF_DOUBLE;

	int num_proposed_moves = 0;
	double accept_rate = 0.0;

	//==== State attributes ====//
	bool _is_start_from_zero_ = false;

	//===== propose : functional methods =====//
	double (Base_Proposer::*_propose_method_)(sys_rnd_gen &, double) = &Base_Proposer::_propose_undefined_;

	// alternative proposer methods
	double _propose_undefined_(sys_rnd_gen &prng, double start)
	{
		basics::raise_error("Attempt to call a method of an undefined proposer.");
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

private:
	//===== State Attributes =====//
	vector<int> _args_to_adapt_ = {};
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  proposer_switchboard  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Base_Proposer *proposer_switchboard(
	string rule_name = UNDEFINED,
	dictargs rule_args = {},
	double rule_from = INF_DOUBLE);

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
		dictargs rule_args = {{"min", (double)-0.5}, {"max", (double)0.5}},
		double rule_from = INF_DOUBLE) : Base_Proposer(new Uniform_Distribution(rule_args), rule_from, {"min", "max"})
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
		dictargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double rule_from = INF_DOUBLE) : Base_Proposer(new Normal_Distribution(rule_args), rule_from, {1})
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
		dictargs rule_args = {{"mean", (double)0.0}, {"sd", (double)1.0}},
		double rule_from = INF_DOUBLE) : Base_Proposer(new Lognorm_Distribution(rule_args), rule_from, {1})
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

template <typename dictargs_val_or_ptr>
class Base_Variable
{
protected:
	string name = UNDEFINED;				// name of the variable
	double value = INF_DOUBLE;				// value of the variable
	double buffer = INF_DOUBLE;				// buffer of the variable
	bool infer = false;						// to be inferred or immutable
	Base_Distribution *distr_ptr = nullptr; // prior/model distribution of the parameter
	Base_Proposer *prop_ptr = nullptr;		// proposer for Monte Carlo Markov Chain

public:
	//=====Constructor=====//
	Base_Variable(const string &name = UNDEFINED, double value = INF_DOUBLE,
				  const string &distr_name = UNDEFINED, dictargs_val_or_ptr distr_args = {},
				  const string &rule_name = UNDEFINED, dictargs rule_args = {},
				  dictargs range = {}, double rule_from = INF_DOUBLE)
	{
		//----- Control options -----//
		// check that distribution args have no mismatch
		if (distr_args.find(IO_MISMATCH) != distr_args.end())
		{
			basics::raise_error("Mismatching definitions in arguments of '" + name + "'.");
		}

		// check that rule args have no mismatch
		if (rule_args.find(IO_MISMATCH) != rule_args.end())
		{
			basics::raise_error("Mismatching definitions in rule arguments of '" + name + "'.");
		}

		//----- Assign options -----//
		// set fundamental necessary values
		set_name(name);
		set_value(value);

		// define distribution_ptr
		distr_ptr = distribution_ptr(distr_name, distr_args);

		// assign its range of validity
		distr_ptr->set_range(range);

		// Assing proposer rule with associated state 'infer'.
		if (basics::is_na(rule_name))
		{
			// immutable variable if no rule is defined
			infer = false;
			_propose_method_ = &Base_Variable::_propose_error_;
		}
		else
		{
			infer = true;

			// define proposer rule with functional programming
			if (rule_name == "gibbs")
			{
				prop_ptr = new Prior_Proposer(distr_ptr);
				_propose_method_ = &Base_Variable::_propose_gibbs_sampl_;
			}
			else
			{
				/*
				// FIXME : this warning doesn't work in a template because global warning are forbidden in headers
				// Waiting for throwing exception structure
				if ((rule_from != INF_DOUBLE) && (!valid(rule_from)))
				{
					std::ostringstream oss;
					oss << "Initial value " << rule_from << " outside of value range"
						<< " {min:" << distr_ptr->get_range_min() << ", max:" << distr_ptr->get_range_max() << "}."
						<< " If it was intentional, ignore this message.";
					screen::raise_warning(oss.str(), Invalid_Proposer_Range_Warning);
				}
				*/

				// assigne the pointer to the proposer
				prop_ptr = proposer_switchboard(rule_name, rule_args, rule_from);
				_propose_method_ = is_from_zero() ? &Base_Variable::_propose_from_zero_ : &Base_Variable::_propose_from_previous_;
			}
		}
	}
	virtual ~Base_Variable() {};

	//=====Getter=====//
	string get_name() const { return name; }
	double get_value() const { return value; }
	const double *get_pointer() const { return &value; }
	bool get_infer() const { return infer; }
	bool is_from_zero() const { return (prop_ptr->is_rule_start_from_zero)(); }
	double get_from_value() const { return (prop_ptr->get_from_value)(); }

	//=====Setter=====//
	void set_name(const string &x) { name = x; }
	void set_value(double x)
	{
		buffer = value;
		value = x;
	}
	void restore_value() { value = buffer; }

	//===== distr_ptr :  =====//
	void distr_rand(sys_rnd_gen &prng) { set_value(distr_ptr->rand(prng)); }
	double distr_logpdf() const { return distr_ptr->logpdf(value); }
	double distr_logratio() const { return distr_ptr->logratio(buffer, value); }
	bool valid() const { return distr_ptr->valid(value); }

	//===== prop_ptr : propose methods =====//
	void propose(sys_rnd_gen &prng) { return set_value((this->*_propose_method_)(prng)); }
	double propose_logratio() const { return (prop_ptr->logratio)(buffer, value); }

	//===== prop_ptr : adaptive methods =====//
	double get_accept_rate() const { return prop_ptr->get_accept_rate(); }
	void update_accept_rate(bool is_accepted) { prop_ptr->update_accept_rate(is_accepted); }
	void reset_accept_rate() { prop_ptr->reset_accept_rate(); }
	void adapt(double ratio) { (prop_ptr->adapt)(ratio); }

	//===== Show =====//
	string show_distr_name() const { return distr_ptr->get_name(); }
	string show_distr_arg_names() const { return distr_ptr->show_arg_names(); }
	string show_distr_arg_values() const { return distr_ptr->show_arg_values(); }

	string show_prop_name() const { return prop_ptr->get_name(); }
	string show_prop_arg_names() const { return prop_ptr->show_arg_names(); }
	string show_prop_arg_values() const { return prop_ptr->show_arg_values(); }

	double get_range_min() const { return distr_ptr->get_range_min(); }
	double get_range_max() const { return distr_ptr->get_range_max(); }

	string show() const
	{
		using namespace basics;
		std::ostringstream oss;
		oss << get_value() << C_DELIM << truth(infer);
		if (infer)
		{
			oss << C_DELIM << show_distr_name() << C_DELIM << show_distr_arg_names() << C_DELIM << show_distr_arg_values()
				<< C_DELIM << show_prop_name() << C_DELIM << show_prop_arg_names() << C_DELIM << show_prop_arg_values()
				<< C_DELIM << get_range_min() << C_DELIM << get_range_max();
		}
		return oss.str();
	}

private:
	//=====Propose : functional methods=====//
	double (Base_Variable::*_propose_method_)(sys_rnd_gen &);
	double _propose_error_(sys_rnd_gen &prng)
	{
		UNUSED(prng);
		basics::raise_error("Call to proposer of immutable variable '" + get_name() + "'.");
		return 0.;
	}
	double _propose_gibbs_sampl_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, 0); }
	double _propose_from_previous_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, get_value()); }
	double _propose_from_zero_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, get_from_value()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  Parameter : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Parameter : public Base_Variable<dictargs>
{
public:
	//===== Constructor =====//
	// NOTE using rule_from = INF_DOUBLE, i. e. default move _propose_from_previous_
	Parameter(const string &param_name = UNDEFINED, double new_val = INF_DOUBLE,
			  const string &prior_name = UNDEFINED, dictargs prior_args = {},
			  const string &rule_name = UNDEFINED, dictargs rule_args = {},
			  dictargs range = {}) : Base_Variable(param_name, new_val, prior_name, prior_args,
												   rule_name, rule_args, range) {};
	
	//===== Aliases for Readability =====//
	void prior_rand(sys_rnd_gen &prng) { distr_rand(prng); }
	double prior_log_pdf() { return distr_logpdf(); } // FIXME : is this ever used?
	double prior_logratio() { return distr_logratio(); }

	//===== Diagnostics =====//
	ulong dgn_count = 0;
	vector<double> dgn_value = {};
	void init_diagnostics(ulong sample_size);
	void update_diagnostics();
	string show_diagnostics() const;
	vector<double> get_dgn_value() const { return dgn_value; }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Augmented : class  //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Augmented : public Base_Variable<dictargs_ptr>
{
public:
	//===== Constructor =====//
	Augmented(const string &augm_name = UNDEFINED, double value = INF_DOUBLE,
			  const string &model_name = UNDEFINED, dictargs_ptr model_args_ptr = {},
			  const string &rule_name = UNDEFINED, dictargs rule_args = {}, dictargs range = {},
			  double rule_from = INF_DOUBLE) : Base_Variable(augm_name, value, model_name, model_args_ptr,
															 rule_name, rule_args, range, rule_from) {};
	
	//===== Aliases for Readability =====//
	void model_rand(sys_rnd_gen &prng) { distr_rand(prng); }
	double model_logpdf() const { return distr_logpdf(); }
	double model_logratio() const { return distr_logratio(); }
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
	vector<string> param_file_headers = {"name", "value", "infer",
										 "prior", "args", "argsval",
										 "rule", "rargs", "rargsval", "min", "max"};
	vector<string> param_to_be_inferred = {};
	vector<string> param_to_be_static = {};

	//===== Constructors =====//
	Monte_Carlo_Core() : _uniform_distr_(0., 1.)
	{
#ifdef _OPENMP
		_max_num_available_threads = omp_get_max_threads();
#endif
	}
	virtual ~Monte_Carlo_Core() {};

	//===== Getter =====//
	void disable_assert_param();
	void assert_param(const string &parID) const;
	double get_param(const string &parID) const { return (this->*_get_param_method_)(parID); }
	const double *get_param_pointer(const string &parID) const { return (this->*_get_param_pointer_method_)(parID); }
	const double *get_reference() const { return &REFERENCE_VALUE; }
	bool get_infer(const string &parID) const { return params.at(parID).get_infer(); }
	double get_accept_rate(const string &parID) const { return params.at(parID).get_accept_rate(); }

	//===== Setter =====//
	void init_rand_params(bool is_rand_params);
	void reset_accept_rate(const string &parID) { params.at(parID).reset_accept_rate(); }
	void assign_param(const string &parID, const Parameter &new_param);
	void drop_all_params();
	const double *assign_scalar_pointer(double new_val);

	//===== Parallelization =====//
	void set_num_threads(unsigned n = 1);
	int get_num_threads() const { return setup_num_threads; };

	//===== Metropolis-Hastings moves =====//
	bool propose_valid_param_move(const string &parID);
	double get_move_logratio() const { return buffer_logratio; }
	void accept_param_move(const string &parID);
	void refuse_param_move(const string &parID);
	void adapt(const string &parID, double ratio) { params.at(parID).adapt(ratio); }
	void init_param_diagnostics(ulong sample_size);
	void update_param_diagnostics();

	//===== Pseudo-Random Number Generator =====//
	void set_prng_seed(int seed);
	sys_rnd_gen &get_prng() { return (this->*_prng_method_)(); }
	double uniform() { return _uniform_distr_(get_prng()); }

	//===== Show =====//
	string show_param_distr_name(const string &parID) const { return params.at(parID).show_distr_name(); }
	string show_param_distr_arg_names(const string &parID) const { return params.at(parID).show_distr_arg_names(); }
	string show_param_distr_arg_values(const string &parID) const { return params.at(parID).show_distr_arg_values(); }
	string show_param_prop_name(const string &parID) const { return params.at(parID).show_prop_name(); }
	string show_param_prop_arg_names(const string &parID) const { return params.at(parID).show_prop_arg_names(); }
	string show_param_prop_arg_values(const string &parID) const { return params.at(parID).show_prop_arg_values(); }
	double get_param_range_min(const string &parID) const { return params.at(parID).get_range_min(); }
	double get_param_range_max(const string &parID) const { return params.at(parID).get_range_max(); }

	string show_param(const string &parID) const;
	string show_param_diagnostics(const string &parID) const;
	vector<double> get_dgn_value(const string &parID) const { return params.at(parID).get_dgn_value(); }

private:
	//===== Internal Attributes =====//
	std::unordered_map<string, Parameter> params = {};	   // parameters dictionary
	std::map<double, Parameter> _scalar_type_params_ = {}; // scalar types dictionary // FIXME : temporary

	double buffer_logratio = 0.0;							// buffer logratio (prior+proposal) of the current parameter
	std::uniform_real_distribution<double> _uniform_distr_; // uniform distribution
	unsigned setup_num_threads = 1;							// number of OMP threads

	//===== Getter (auxiliary methods) =====//
	double (Monte_Carlo_Core::*_get_param_method_)(const string &) const = &Monte_Carlo_Core::_get_param_assert_;
	double _get_param_assert_(const string &parID) const;
	double _get_param_default_(const string &parID) const;

	const double *(Monte_Carlo_Core::*_get_param_pointer_method_)(const string &) const = &Monte_Carlo_Core::_get_param_pointer_assert_;
	const double *_get_param_pointer_assert_(const string &parID) const;
	const double *_get_param_pointer_default_(const string &parID) const;

	//===== Pseudo Random Number Generator (auxiliary methods) =====//
	// FIXME : this should change if Rcpp according to developer guides
	sys_rnd_gen prng;						// Mersenne Twister 19937 generator (32 or 64 bit)
	vector<sys_rnd_gen> prng_parallel = {}; // parallel generators (for multi-threading)
	sys_rnd_gen &(Monte_Carlo_Core::*_prng_method_)();
	sys_rnd_gen &_prng_unique_thread_() { return prng; }
#ifdef _OPENMP
	unsigned _max_num_available_threads = 1;
	sys_rnd_gen &_prng_parallel_threads_() { return prng_parallel[omp_get_thread_num()]; }
#endif
};

#endif // VARIA__H