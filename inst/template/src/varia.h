/*
 * varia.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef VARIA__H
#define VARIA__H

#include <memory>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>

#include "tools.h"
#include "function.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* A parent class to all distributions. */
class Base_Distribution : public Base_Function
{
public:
	using Base_Function::Base_Function;

	double pdf(double x) const { return (this->*_func_pdf_)(x); }
	double cdf(double x) const { return (this->*_func_cdf_)(x); }
	double logpdf(double x) const { return (this->*_func_logpdf_)(x); }
	double logratio(double xold, double xnew) const { return (this->*_func_logratio_)(xold, xnew); }
	double rand(sys_rnd_gen &prng) { return (this->*_func_rand_)(prng); }

protected:
	//----- Distribution specific virtual methods -----//
	virtual double _pdf_(double x) const;
	virtual double _cdf_(double x) const;
	virtual double _logpdf_(double x) const;
	virtual double _logratio_(double xold, double xnew) const;
	virtual double _rand_(sys_rnd_gen &prng);

	//----- Adjust methods for truncation and/or discretizat -----//
	void _assign_reshaping_(const dictargs &range, const std::string &discr);
	void _check_rand_attempts_(int &num_attempts) const;

	//----- Functional methods for public methods -----//
	double (Base_Distribution::*_func_pdf_)(double) const = &Base_Distribution::_default_pdf_;
	double (Base_Distribution::*_func_cdf_)(double) const = &Base_Distribution::_default_cdf_;
	double (Base_Distribution::*_func_logpdf_)(double) const = &Base_Distribution::_default_logpdf_;
	double (Base_Distribution::*_func_logratio_)(double, double) const = &Base_Distribution::_default_logratio_;
	double (Base_Distribution::*_func_rand_)(sys_rnd_gen &) = &Base_Distribution::_default_rand_;

	//----- Functional methods for truncation -----//
	double (Base_Distribution::*_func_trunc_)(double) const = &Base_Distribution::_trunc_dummy_;
	double (Base_Distribution::*_func_logtrunc_)(double) const = &Base_Distribution::_logtrunc_dummy_;
	double (Base_Distribution::*_func_shift_)(double) const = &Base_Distribution::_shift_dummy_;

private:
	double _undefined_call_error_(const std::string &txt = UNDEFINED, double x = 0, double y = 0) const;

	//----- Deafult functional methods -----//
	double _default_pdf_(double x) const;
	double _default_cdf_(double x) const;
	double _default_logpdf_(double x) const;
	double _default_logratio_(double xold, double xnew) const;
	double _default_rand_(sys_rnd_gen &prng);

	//----- Functional methods for the center adjusted discretization -----//
	double _discr_center_pdf_(double x) const;
	double _discr_center_cdf_(double x) const;
	double _discr_center_logpdf_(double x) const;
	double _discr_center_logratio_(double xold, double xnew) const;
	double _discr_center_rand_(sys_rnd_gen &prng);

	//----- Functional methods for the left adjusted discretization -----//
	double _discr_left_pdf_(double x) const;
	double _discr_left_cdf_(double x) const;
	double _discr_left_logpdf_(double x) const;
	double _discr_left_logratio_(double xold, double xnew) const;
	double _discr_left_rand_(sys_rnd_gen &prng);

	//----- Functional methods for the truncation and/or discretization -----//
	double _trunc_dummy_(double x) const;
	double _trunc_min_(double x) const;
	double _trunc_max_(double x) const;
	double _trunc_both_sides_(double x) const;

	//----- Functional methods for the (log) truncation -----//
	double _logtrunc_dummy_(double x) const;
	double _logtrunc_min_(double x) const;
	double _logtrunc_max_(double x) const;
	double _logtrunc_both_sides_(double x) const;

	//----- Functional methods to shift cumulative -----//
	double _shift_dummy_(double x) const;
	double _shift_min_(double x) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  distribution switchboard  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* overloaded distribution with static values. */
std::shared_ptr<Base_Distribution> distribution_ptr(const std::string &name = UNDEFINED,
													  const dictargs &args = {},
													  const dictargs &range = {});

/* overloaded distribution with pointer to values. */
std::shared_ptr<Base_Distribution> distribution_ptr(const std::string &name = UNDEFINED,
													  const dictargs_ptr &args = {},
													  const dictargs &range = {});

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
	Base_Proposer(std::shared_ptr<Base_Distribution>, double rule_start_value = INF_DOUBLE,
					std::vector<int> args_to_adapt = {});
	virtual ~Base_Proposer() {};

	//===== Getter =====//
	const std::string &get_name() const { return distr_ptr->get_name(); }
	double get_arg(int argID) const { return (distr_ptr->get_arg)(argID); }
	double get_from_value() const { return from_value; }
	bool is_rule_start_from_zero() const { return _is_start_from_zero_; }

	//===== Setter =====//
	void set_rule_start_value(double x) { from_value = x; }

	double rand(sys_rnd_gen &prng) { return (distr_ptr->rand)(prng); }

	virtual double propose(sys_rnd_gen &prng, double start = 0.0) { return (this->*_func_propose_)(prng, start); }
	virtual double logratio(double xold, double xnew) const { return (this->*_logratio_method_)(xold, xnew); }

	//=====Finite Adaptive=====//
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

	//===== Show =====//
	const std::string show_arg_names() const { return distr_ptr->show_arg_names(); }
	const std::string show_arg_values() const { return distr_ptr->show_arg_values(); }

protected:
	//===== Attributes =====//
	std::shared_ptr<Base_Distribution> distr_ptr;
	double from_value = INF_DOUBLE;

	int num_proposed_moves = 0;
	double accept_rate = 0.0;

	//==== State attributes ====//
	bool _is_start_from_zero_ = false;

	//===== propose : functional methods =====//
	double (Base_Proposer::*_func_propose_)(sys_rnd_gen &, double) = &Base_Proposer::_propose_undefined_;

	// alternative proposer methods
	double _propose_undefined_(sys_rnd_gen &prng, double start)
	{
		tools::error("Attempt to call a method of an undefined proposer.");
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

	//===== State Attributes =====//
	std::vector<int> _args_to_adapt_ = {};
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

class Base_Variable
{
protected:
	std::string name = UNDEFINED;							  // name of the variable
	double value = INF_DOUBLE;								  // value of the variable
	double buffer = INF_DOUBLE;								  // buffer of the variable
	bool infer = false;										  // to be inferred or immutable
	std::shared_ptr<Base_Distribution> distr_ptr = nullptr; // prior/model distribution of the parameter
	std::shared_ptr<Base_Proposer> prop_ptr = nullptr;	  // proposer for Monte Carlo Markov Chain

public:
	//=====Constructor=====//
	Base_Variable(const std::string &name_ = UNDEFINED,
				  double value_ = std::numeric_limits<double>::quiet_NaN()) : name(name_), value(value_) {}
	virtual ~Base_Variable() {};

	// define a pointer to distribution
	template <typename dictargs_val_or_ptr>
	void assign_distribution(const std::string &distr_name = UNDEFINED,
							 dictargs_val_or_ptr distr_args = {},
							 dictargs range = {})
	{
		// check that distribution args have no mismatch
		if (distr_args.find(IO_MISMATCH) != distr_args.end())
		{
			tools::error("Mismatching definitions in arguments of '" + name + "'.");
		}
		distr_ptr = distribution_ptr(distr_name, distr_args, range);
	}

	// Assing a proposer rule with associated state 'infer'.
	void assign_proposer(const std::string &rule_name = UNDEFINED,
						 dictargs rule_args = {},
						 double rule_from = INF_DOUBLE);

	//===== Getter =====//const
	const std::string &get_name() const { return name; }
	double get_value() const { return value; }
	const double *get_pointer() const { return &value; }
	bool get_infer() const { return infer; }
	bool is_from_zero() const { return (prop_ptr->is_rule_start_from_zero)(); }
	double get_from_value() const { return (prop_ptr->get_from_value)(); }

	//===== Setter =====//
	void set_name(const std::string &x) { name = x; }
	void set_value(double x)
	{
		buffer = value;
		value = x;
	}
	void restore_value() { value = buffer; }

	//===== Metropolis-Hastings =====//

	double get_move_logratio() const { return (this->*_func_move_logratio_)(); }
	bool valid_move() const { return (this->*_func_valid_move_)(); }

	//----- distr_ptr : prior/model methods -----//
	void rand(sys_rnd_gen &prng) { set_value(distr_ptr->rand(prng)); }

	double distr_logpdf() const { return distr_ptr->logpdf(value); }
	double distr_logratio() const { return distr_ptr->logratio(buffer, value); }

	//----- prop_ptr : propose methods -----//
	void propose(sys_rnd_gen &prng) { return set_value((this->*_func_propose_)(prng)); }
	double propose_logratio() const { return (prop_ptr->logratio)(buffer, value); }

	//----- prop_ptr : adaptive methods -----//
	double get_accept_rate() const { return prop_ptr->get_accept_rate(); }
	void update_accept_rate(bool is_accepted) { prop_ptr->update_accept_rate(is_accepted); }
	void reset_accept_rate() { prop_ptr->reset_accept_rate(); }
	void adapt_rule(double ratio) { (prop_ptr->adapt)(ratio); }

	//===== Show =====//
	const std::string show_distr_name() const { return distr_ptr->get_name(); }
	const std::string show_distr_arg_names() const { return distr_ptr->show_arg_names(); }
	const std::string show_distr_arg_values() const { return distr_ptr->show_arg_values(); }

	const std::string show_prop_name() const { return prop_ptr->get_name(); }
	const std::string show_prop_arg_names() const { return prop_ptr->show_arg_names(); }
	const std::string show_prop_arg_values() const { return prop_ptr->show_arg_values(); }

	const std::string show_range_min() const { return distr_ptr->show_range_min(); }
	const std::string show_range_max() const { return distr_ptr->show_range_max(); }

	const std::string show() const
	{
		using namespace tools;
		std::ostringstream oss;
		oss << get_value() << C_DELIM << truth(infer);
		if (infer)
		{
			oss << C_DELIM << show_distr_name() << C_DELIM << show_distr_arg_names() << C_DELIM << show_distr_arg_values()
				<< C_DELIM << show_prop_name() << C_DELIM << show_prop_arg_names() << C_DELIM << show_prop_arg_values()
				<< C_DELIM << (tools::is_na(show_range_min()) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(show_range_min()))
				<< C_DELIM << (tools::is_na(show_range_max()) ? std::numeric_limits<double>::quiet_NaN() : tools::stod(show_range_max()));
		}
		return oss.str();
	}

private:
	//===== Valid move : functional methods =====//
	bool (Base_Variable::*_func_valid_move_)() const = &Base_Variable::_valid_move_default_;
	bool _valid_move_default_() const { return distr_ptr->valid(value); }
	bool _valid_move_gibbs_() const { return true; }

	//===== Move logratio : functional methods =====//
	double (Base_Variable::*_func_move_logratio_)() const = &Base_Variable::_move_logratio_default_;
	double _move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
	double _move_logratio_gibbs_() const { return 0.0; }

	//===== Propose : functional methods =====//
	double (Base_Variable::*_func_propose_)(sys_rnd_gen &);
	double _propose_error_(sys_rnd_gen &prng)
	{
		UNUSED(prng);
		tools::error("Call to proposer of immutable variable '" + get_name() + "'.");
		return 0.;
	}
	double _propose_gibbs_sample_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, 0); }
	double _propose_from_previous_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, get_value()); }
	double _propose_from_zero_(sys_rnd_gen &prng) { return (prop_ptr->propose)(prng, get_from_value()); }
};

#endif // VARIA__H