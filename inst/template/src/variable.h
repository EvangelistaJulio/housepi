/*
 * variable.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__VARIABLE__H
#define HOUSEPI__VARIABLE__H

#include <memory>
#include <string>
#include <vector>

#include "tools.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Variable : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Variable
{
public:
	//===== CONSTRUCTOR =====//
	Variable(const std::string &name = NULL_STRING,
				  double value = std::numeric_limits<double>::quiet_NaN(),
				  bool is_discrete = false);
	Variable(Variable &&other) noexcept;			  // Move constructor
	Variable &operator=(Variable &&other) noexcept; // Move assignment
	~Variable();

	//===== GETTER =====//
	double get_value() const;
	const double *get_pointer() const;
	double get_from_value() const;
	bool is_constant() const;
	bool is_discrete() const;
	bool is_from_zero() const;

	//===== SETTER =====//
	void init(const std::string &name, double value, bool is_discrete);
	void user_reset();
	void set_value(double x);
	void restore_value();

	//----- Pointer to distribution -----//
	void assign_distribution(const std::string &distr_name = NULL_STRING);

	void assign_distr_scalar_arg(const std::string &arg_name, double x);
	void assign_distr_pointer_arg(const std::string &arg_name, const double *y);

	void assign_distr_scalar_range_min(double x);
	void assign_distr_pointer_range_min(const double *y);

	void assign_distr_scalar_range_max(double x);
	void assign_distr_pointer_range_max(const double *y);

	void init_distr_reshaping();

	//----- Pointer to proposer -----//
	void assign_proposer(const std::string &rule_name = NULL_STRING,
						 double rule_scale = NULL_DOUBLE,
						 double rule_from = NULL_DOUBLE);

	//===== METROPOLIS-HASTINGS =====//

	double get_move_logratio() const;
	bool valid_move() const;

	//----- distr_ptr : prior/model methods -----//
	void rand(rnd_gen &prng);

	double distr_logpdf() const;
	double distr_logratio() const;

	//----- prop_ptr : propose methods -----//
	void propose(rnd_gen &prng);
	double propose_logratio() const;

	//----- prop_ptr : adaptive methods -----//
	double get_accept_rate() const;
	void update_accept_rate(bool is_accepted);
	void reset_accept_rate();
	void adapt_rule(double ratio);

	//===== SHOW =====//
	const std::string show_name() const;

	const std::string show_distr_name() const;
	const std::string show_distr_arg_names() const;
	const std::string show_distr_arg_values() const;

	const std::string show_min() const;
	const std::string show_max() const;

	const std::string show_prop_name() const;
	const std::string show_prop_arg_names() const;
	const std::string show_prop_arg_values() const;
	const std::string show_prop_scale() const;

private:
	//===== Pointer to Implementation =====//
	class Impl;
	std::unique_ptr<Impl> pImpl;
};

#endif // HOUSEPI__VARIABLE__H
