/*
 * variable.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__VARIABLE__H
#define HOUSEPI__VARIABLE__H

#include <memory>
#include <random>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>

#include "tools.h"
#include "distribution.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base_Proposer : class  /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Base_Proposer
{
public:
	//===== CONSTRUCTOR =====//
	Base_Proposer(std::shared_ptr<Base_Distribution>,
				  double rule_start_value = INF_DOUBLE,
				  std::vector<ushort> args_to_adapt = {});
	virtual ~Base_Proposer() {};

	//===== GETTER =====//
	const std::string &get_name() const;
	double get_arg(ushort argID) const;
	double get_from_value() const;
	bool is_rule_start_from_zero() const;

	//===== SETTER =====//
	void set_rule_start_value(double x);
	double rand(sys_rnd_gen &prng);
	virtual double propose(sys_rnd_gen &prng, double start = 0.0);
	virtual double logratio(double xold, double xnew) const;

	//===== ADAPTATION =====//
	double get_accept_rate() const;
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt(double ratio);

	//===== SHOW =====//
	const std::string show_arg_names() const;
	const std::string show_arg_values() const;

protected:
	std::shared_ptr<Base_Distribution> distr_ptr;
	double from_value = INF_DOUBLE;

	ulong n_proposed_moves = 0;
	double accept_rate = 0.0;

	bool _is_start_from_zero_ = false;
	std::vector<ushort> _args_to_adapt_ = {};

	//===== propose : functional methods =====//
	double (Base_Proposer::*_func_propose_)(sys_rnd_gen &, double) = &Base_Proposer::_propose_undefined_;
	double _propose_undefined_(sys_rnd_gen &prng, double start);
	double _propose_sum_to_previous_(sys_rnd_gen &prng, double start);
	double _propose_times_previous_(sys_rnd_gen &prng, double start);
	double _propose_draw_rand_(sys_rnd_gen &prng, double start);

	//===== logratio : functional methods =====//
	double (Base_Proposer::*_logratio_method_)(double, double) const = &Base_Proposer::_logratio_undefined_;
	double _logratio_undefined_(double xold, double xnew) const;
	double _logratio_pdf_(double xold, double xnew) const;
	double _logratio_sum_to_zero_(double xold, double xnew) const;
	double _logratio_times_zero_(double xold, double xnew) const;
	double _logratio_symm_(double xold, double xnew) const;
	double _logratio_lognorm_symm_(double xold, double xnew) const;
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
	std::string name = UNDEFINED; // name of the variable
	double value = INF_DOUBLE;	  // value of the variable
	double buffer = INF_DOUBLE;	  // buffer of the variable
	bool is_discr = false;        // wheter the variable is discrete or not
	bool infer = false;			  // to be inferred or immutable

	std::shared_ptr<Base_Distribution> distr_ptr = nullptr; // prior/model distribution of the variable
	std::shared_ptr<Base_Proposer> prop_ptr = nullptr;		// proposer for Monte Carlo

public:
	//===== CONSTRUCTOR =====//
	Base_Variable(const std::string &name = UNDEFINED,
				  double value = std::numeric_limits<double>::quiet_NaN(),
				  bool is_discr = false);
	virtual ~Base_Variable() {};

	//===== GETTER =====//
	const std::string &get_name() const;
	double get_value() const;
	const double *get_pointer() const;
	bool get_infer() const;
	bool is_from_zero() const;
	double get_from_value() const;

	//===== SETTER =====//
	void set_name(const std::string &x);
	void set_value(double x);
	void restore_value();

	// define a pointer to distribution
	void assign_distribution(const std::string &distr_name = UNDEFINED,
							 dictargs distr_args = {},
							 dictargs range = {});
	void assign_distribution(const std::string &distr_name = UNDEFINED,
							 dictargs_ptr distr_args = {},
							 dictargs range = {});

	// Assing a proposer rule with associated state 'infer'.
	void assign_proposer(const std::string &rule_name = UNDEFINED,
						 double rule_args = INF_DOUBLE,
						 double rule_from = INF_DOUBLE);

	//===== METROPOLIS-HASTINGS =====//

	double get_move_logratio() const;
	bool valid_move() const;

	//----- distr_ptr : prior/model methods -----//
	void rand(sys_rnd_gen &prng);

	double distr_logpdf() const;
	double distr_logratio() const;

	//----- prop_ptr : propose methods -----//
	void propose(sys_rnd_gen &prng);
	double propose_logratio() const;

	//----- prop_ptr : adaptive methods -----//
	double get_accept_rate() const;
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt_rule(double ratio);

	//===== SHOW =====//
	const std::string show_distr_name() const;
	const std::string show_distr_arg_names() const;
	const std::string show_distr_arg_values() const;

	const std::string show_prop_name() const;
	const std::string show_prop_arg_names() const;
	const std::string show_prop_arg_values() const;

	const std::string show_range_min() const;
	const std::string show_range_max() const;

	const std::string show() const;

private:
	//===== Valid move : functional methods =====//
	bool (Base_Variable::*_func_valid_move_)() const = &Base_Variable::_valid_move_default_;
	bool _valid_move_default_() const;
	bool _valid_move_gibbs_() const;

	//===== Move logratio : functional methods =====//
	double (Base_Variable::*_func_move_logratio_)() const = &Base_Variable::_move_logratio_default_;
	double _move_logratio_default_() const;
	double _move_logratio_gibbs_() const;

	//===== Propose : functional methods =====//
	double (Base_Variable::*_func_propose_)(sys_rnd_gen &) = &Base_Variable::_propose_error_;
	double _propose_error_(sys_rnd_gen &prng);
	double _propose_gibbs_sample_(sys_rnd_gen &prng);
	double _propose_from_previous_(sys_rnd_gen &prng);
	double _propose_from_zero_(sys_rnd_gen &prng);
};

#endif // HOUSEPI__VARIABLE__H
