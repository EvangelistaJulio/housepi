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

#include "utils.h"

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

    double pdf(double x) const { return _trunc_(_pdf_(x)); }
    double cdf(double x) const { return _trunc_(_shift_(_cdf_(x))); }
    double logpdf(double x) const { return _logtrunc_(_logpdf_(x)); }

    //----- Distribution specific methods -----//
    virtual double logratio(double xold, double xnew) const
    {
        UNUSED(xold);
        UNUSED(xnew);
        return 0.0;
    }
    virtual double rand(sys_rnd_gen &prng)
    {
        UNUSED(prng);
        return _undefined_call_error_("rand");
    }

protected:
    //----- Distribution specific methods -----//
    virtual double _pdf_(double x) const { return _undefined_call_error_("pdf", x); }
    virtual double _logpdf_(double x) const { return _undefined_call_error_("logpdf", x); }
    virtual double _cdf_(double x) const { return _undefined_call_error_("cdf", x); }

    //----- Adjust methods for truncation -----//
    void set_truncation(const dictargs &range);
    void check_rand_attempts(int &num_attempts) const;

private:
    //----- Internal error for calls to undefined -----//
    double _undefined_call_error_(const std::string &txt, double x = 0, double y = 0) const;

    //----- Functional methods for the truncation -----//
    double (Base_Distribution::*_trunc_method_)(double) const = &Base_Distribution::_trunc_dummy_;
    double _trunc_(double x) const { return (this->*_trunc_method_)(x); }
    double _trunc_dummy_(double x) const { return x; }
    double _trunc_min_(double x) const { return x / (1.0 - _cdf_(get_range_min())); }
    double _trunc_max_(double x) const { return x / (_cdf_(get_range_max())); }
    double _trunc_both_sides_(double x) const { return x / (_cdf_(get_range_max()) - _cdf_(get_range_max())); }

    //----- Functional methods for the (log) truncation -----//
    double (Base_Distribution::*_logtrunc_method_)(double) const = &Base_Distribution::_logtrunc_dummy_;
    double _logtrunc_(double x) const { return (this->*_logtrunc_method_)(x); }
    double _logtrunc_dummy_(double x) const { return x; }
    double _logtrunc_min_(double x) const { return x - log(1.0 - _cdf_(get_range_min())); }
    double _logtrunc_max_(double x) const { return x - log(_cdf_(get_range_max())); }
    double _logtrunc_both_sides_(double x) const { return x - log(_cdf_(get_range_max()) - _cdf_(get_range_min())); }

    //----- Functional methods to shift cumulative -----//
    double (Base_Distribution::*_shift_method_)(double) const = &Base_Distribution::_shift_dummy_;
    double _shift_(double x) const { return (this->*_shift_method_)(x); }
    double _shift_dummy_(double x) const { return x; }
    double _shift_min_(double x) const { return x - _cdf_(get_range_min()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  distribution switchboard  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

Base_Distribution *distribution_ptr(const std::string &name = UNDEFINED,
                                    const dictargs &args = {},
                                    const dictargs &range = {});
Base_Distribution *distribution_ptr(const std::string &name = UNDEFINED,
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
	Base_Proposer(Base_Distribution *, double rule_start_value = INF_DOUBLE,
				  std::vector<int> args_to_adapt = {});
	virtual ~Base_Proposer() {};

	//===== Getter =====//
	std::string get_name() const { return distr_ptr->get_name(); }
	double get_arg(int argID) const { return (distr_ptr->get_arg)(argID); }
	double get_from_value() const { return from_value; }
	bool is_rule_start_from_zero() const { return _is_start_from_zero_; }

	//===== Setter =====//
	void set_rule_start_value(double x) { from_value = x; }

	double rand(sys_rnd_gen &prng) { return (distr_ptr->rand)(prng); }

	virtual double propose(sys_rnd_gen &prng, double start = 0.0) { return (this->*_propose_method_)(prng, start); }
	virtual double logratio(double xold, double xnew) const { return (this->*_logratio_method_)(xold, xnew); }

	//=====Finite Adaptive=====//
	double get_accept_rate() const { return accept_rate; };
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

	//===== Show =====//
	std::string show_arg_names() const { return distr_ptr->show_arg_names(); }
	std::string show_arg_values() const { return distr_ptr->show_arg_values(); }

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
	std::string name = UNDEFINED;			// name of the variable
	double value = INF_DOUBLE;				// value of the variable
	double buffer = INF_DOUBLE;				// buffer of the variable
	bool infer = false;						// to be inferred or immutable
	Base_Distribution *distr_ptr = nullptr; // prior/model distribution of the parameter
	Base_Proposer *prop_ptr = nullptr;		// proposer for Monte Carlo Markov Chain

public:
	//=====Constructor=====//
	Base_Variable(const std::string &name_, double value_) : name(name_), value(value_) {}
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

	//===== Getter =====//
	std::string get_name() const { return name; }
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

	double get_move_logratio() const { return (this->*_move_logratio_method_)(); }
	bool valid_move() const { return (this->*_valid_move_method_)(); }

	//----- distr_ptr : prior/model methods -----//
	void rand(sys_rnd_gen &prng) { set_value(distr_ptr->rand(prng)); }

	double distr_logpdf() const { return distr_ptr->logpdf(value); }
	double distr_logratio() const { return distr_ptr->logratio(buffer, value); }

	//----- prop_ptr : propose methods -----//
	void propose(sys_rnd_gen &prng) { return set_value((this->*_propose_method_)(prng)); }
	double propose_logratio() const { return (prop_ptr->logratio)(buffer, value); }

	//----- prop_ptr : adaptive methods -----//
	double get_accept_rate() const { return prop_ptr->get_accept_rate(); }
	void update_accept_rate(bool is_accepted) { prop_ptr->update_accept_rate(is_accepted); }
	void reset_accept_rate() { prop_ptr->reset_accept_rate(); }
	void adapt_rule(double ratio) { (prop_ptr->adapt)(ratio); }

	//===== Show =====//
	std::string show_distr_name() const { return distr_ptr->get_name(); }
	std::string show_distr_arg_names() const { return distr_ptr->show_arg_names(); }
	std::string show_distr_arg_values() const { return distr_ptr->show_arg_values(); }

	std::string show_prop_name() const { return prop_ptr->get_name(); }
	std::string show_prop_arg_names() const { return prop_ptr->show_arg_names(); }
	std::string show_prop_arg_values() const { return prop_ptr->show_arg_values(); }

	std::string show_range_min() const { return distr_ptr->show_range_min(); }
	std::string show_range_max() const { return distr_ptr->show_range_max(); }

	std::string show() const
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
	bool (Base_Variable::*_valid_move_method_)() const = &Base_Variable::_valid_move_default_;
	bool _valid_move_default_() const { return distr_ptr->valid(value); }
	bool _valid_move_gibbs_() const { return true; }

	//===== Move logratio : functional methods =====//
	double (Base_Variable::*_move_logratio_method_)() const = &Base_Variable::_move_logratio_default_;
	double _move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
	double _move_logratio_gibbs_() const { return 0.0; }

	//===== Propose : functional methods =====//
	double (Base_Variable::*_propose_method_)(sys_rnd_gen &);
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