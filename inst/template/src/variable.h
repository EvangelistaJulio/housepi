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
#include "distribution.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Proposer : class  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Base class for Metropolis proposers. */
class Proposer
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//
	
    Proposer() {};
    Proposer(std::shared_ptr<Distribution>,
             double rule_scale = NULL_DOUBLE,
             double rule_start_value = NULL_DOUBLE,
             const std::vector<ushort> &args_to_adapt = {},
             bool is_discrete = false);
    virtual ~Proposer() {};

    //==================//
    //===== GETTER =====//
    //==================//

    const std::string show_name() const;
    double get_arg(ushort argID) const;
    double get_from_value() const;
    double propose(rnd_gen &prng, double start = 0.);
    double logratio(double xold, double xnew) const;
    double rand(rnd_gen &prng);
    bool is_rule_start_from_zero() const;

    //==================//
    //===== SETTER =====//
    //==================//

    virtual void user_reset();
    void set_rule_start_value(double x);

    //======================//
    //===== ADAPTATION =====//
    //======================//

    double get_accept_rate() const;
    void update_accept_rate(bool is_accepted);
    void reset_accept_rate();
    void adapt(double ratio);

    //================//
    //===== SHOW =====//
    //================//

    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_scale() const;

protected:
    std::string name = NULL_STRING;
    std::shared_ptr<Distribution> _distr_ptr_;

    double _user_rule_scale_ = 1.;     // scale of the Metropolis rule
    bool _is_start_from_zero_ = false; // wheter to move from given value
    double from_value = NULL_DOUBLE;   // given value where to start from

    double accept_rate = 0.;    // acceptance rate
    ulong n_proposed_moves = 0; // number of proposed moves

    ushort n_args_to_adapt = 0;        // how many arguments to adapt
    ushort *_args_to_adapt_ = nullptr; // which arguments should be adapted

    //===== propose : functional methods =====//
    double (Proposer::*_func_propose_)(rnd_gen &, double) = &Proposer::_propose_undefined_;
    double _propose_undefined_(rnd_gen &prng, double start);
    double _propose_sum_to_previous_(rnd_gen &prng, double start);
    double _propose_times_previous_(rnd_gen &prng, double start);
    double _propose_draw_rand_(rnd_gen &prng, double start);

    //===== logratio : functional methods =====//
    double (Proposer::*_func_logratio_)(double, double) const = &Proposer::_logratio_undefined_;
    double _logratio_undefined_(double xold, double xnew) const;
    double _logratio_pdf_(double xold, double xnew) const;
    double _logratio_sum_to_zero_(double xold, double xnew) const;
    double _logratio_times_zero_(double xold, double xnew) const;
    double _logratio_symm_(double xold, double xnew) const;
    double _logratio_lognorm_symm_(double xold, double xnew) const;
};

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
	~Variable();

	//===== GETTER =====//
	double get() const;
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
    std::shared_ptr<Distribution> _distr_ptr_; // prior\model distribution of the variable
    std::shared_ptr<Proposer> _prop_ptr_;      // Monte Carlo proposer

    std::string name = NULL_STRING; // name of the variable
    double value = NULL_DOUBLE;     // value of the variable
    double buffer = NULL_DOUBLE;    // buffer of the variable

    double _user_value_ = NULL_DOUBLE; // initial user defined value
    bool _is_discrete_ = false;        // wheter discrete or not
    bool _is_constant_ = true;         // wheter constant or not

    //===============================//
    //===== METROPOLIS-HASTINGS =====//
    //===============================//

    //----- Valid move : functional methods -----//
    bool (Variable::*_func_valid_move_)() const = &Variable::_valid_move_default_;
    bool _valid_move_default_() const;
    bool _valid_move_gibbs_() const;

    //----- Move logratio : functional methods -----//
    double (Variable::*_func_move_logratio_)() const = &Variable::_move_logratio_default_;
    double _move_logratio_default_() const;
    double _move_logratio_gibbs_() const;

    //----- Propose : functional methods -----//
    double (Variable::*_func_propose_)(rnd_gen &) = &Variable::_propose_error_;
    double _propose_error_(rnd_gen &prng);
    double _propose_gibbs_sample_(rnd_gen &prng);
    double _propose_from_previous_(rnd_gen &prng);
    double _propose_from_zero_(rnd_gen &prng);
};

#endif // HOUSEPI__VARIABLE__H
