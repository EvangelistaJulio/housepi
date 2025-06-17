/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * variable.h
 *
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
             double rule_scale = HOUSEPI_NA_DOUBLE,
             double rule_start_value = HOUSEPI_NA_DOUBLE,
             const std::vector<ushort> &args_to_adapt = {},
             bool is_discrete = false);
    Proposer(const std::string &distr_name,
             double rule_scale = HOUSEPI_NA_DOUBLE,
             double rule_start_value = HOUSEPI_NA_DOUBLE,
             const std::vector<ushort> &args_to_adapt = {},
             bool is_discrete = false);
    virtual ~Proposer() {};

    //==================//
    //===== GETTER =====//
    //==================//

    const std::string show_name() const { return this->_proposer_name; }
    double get_arg(ushort argID) const { return (_distr_ptr->get_arg)(argID); }
    double get_from_value() const { return _from_value; }
    double propose(rnd_gen &prng, double start) { return (this->*_do_propose)(prng, start); }
    double logratio(double xold, double xnew) const { return (this->*_do_logratio)(xold, xnew); }
    double rand(rnd_gen &prng) { return (_distr_ptr->rand)(prng); }
    bool is_rule_start_from_zero() const { return _is_start_from_zero; }

    //==================//
    //===== SETTER =====//
    //==================//

    virtual void user_reset() {}
    void set_rule_start_value(double x) { _from_value = x; }

    //======================//
    //===== ADAPTATION =====//
    //======================//

    double get_accept_rate() const { return _accept_rate; }
    void update_accept_rate(bool is_accepted)
    {
        _accept_rate = (double)(_accept_rate * _n_proposed_moves + is_accepted) / (_n_proposed_moves + 1);
        ++_n_proposed_moves;
    }
    void reset_accept_rate()
    {
        _accept_rate = 0.;
        _n_proposed_moves = 0;
    }
    void adapt(double ratio);

    //================//
    //===== SHOW =====//
    //================//

    const std::string show_arg_names() const { return _distr_ptr->show_arg_names(); }
    const std::string show_arg_values() const { return _distr_ptr->show_arg_values(); }
    const std::string show_scale() const;

protected:
    std::string _proposer_name = HOUSEPI_NA_STRING;
    std::shared_ptr<Distribution> _distr_ptr;

    double _user_rule_scale = 1.;           // scale of the Metropolis rule
    bool _is_start_from_zero = false;       // wheter to move from given value
    double _from_value = HOUSEPI_NA_DOUBLE; // given value where to start from

    double _accept_rate = 0.;    // acceptance rate
    ulong _n_proposed_moves = 0; // number of proposed moves

    ushort _n_args_to_adapt = 0;           // how many arguments to adapt
    ushort *_args_to_adapt_list = nullptr; // which arguments should be adapted

    //===== propose : functional methods =====//
    double (Proposer::*_do_propose)(rnd_gen &, double) = &Proposer::__propose_undefined;
    double __propose_undefined(rnd_gen &, double) { throw std::invalid_argument(
        "Called method of an undefined proposer."); }
    double __propose_sum_to_previous(rnd_gen &prng, double start) { return start + rand(prng); }
    double __propose_times_previous(rnd_gen &prng, double start) { return start * rand(prng); }
    double __propose_draw_rand(rnd_gen &prng, double) { return rand(prng); }

    //===== logratio : functional methods =====//
    double (Proposer::*_do_logratio)(double, double) const = &Proposer::__logratio_undefined;
    double __logratio_undefined(double, double) const { return 0.; }
    double __logratio_pdf(double xold, double xnew) const { return _distr_ptr->logratio(xold, xnew); }
    double __logratio_sum_to_zero(double xold, double xnew) const { return __logratio_pdf(xold - _from_value, xnew - _from_value); }
    double __logratio_times_zero(double xold, double xnew) const { return __logratio_pdf(xold / _from_value, xnew / _from_value); }
    double __logratio_symm(double, double) const { return 0.; }
    double __logratio_lognorm_symm(double xold, double xnew) const { return std::log(xold / xnew); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Variable : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Variable
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Variable(const std::string &name = HOUSEPI_NA_STRING,
             double value = HOUSEPI_NA_DOUBLE,
             bool is_discrete = false,
             bool is_constant = true,
             const std::string &func_name = HOUSEPI_NA_STRING)
        : _func_ptr(std::make_shared<Function>()),
          _distr_ptr(std::make_shared<Distribution>()),
          _prop_ptr(std::make_shared<Proposer>())
    {
        init(name, value, is_discrete, is_constant, func_name);
    }
    ~Variable() = default;

    double operator()(double x) const { return (this->*_do_operator_return)(x); }

    //==================//
    //===== GETTER =====//
    //==================//

    double get() const { return (this->*_do_get)(); }
    const double *get_pointer() const { return (this->*_do_get_pointer)(); }
    std::shared_ptr<Function> get_function() const { return _func_ptr; }

    bool is_function() const { return _do_get != &Variable::__get_default; }
    bool is_constant() const { return _is_constant; }
    bool is_discrete() const { return _distr_ptr->is_discrete(); }
    bool is_from_zero() const { return (_prop_ptr->is_rule_start_from_zero)(); }

    double get_from_value() const { return (_prop_ptr->get_from_value)(); }

    //==================//
    //===== SETTER =====//
    //==================//

    void init(const std::string &name, double value, bool is_discrete,
              bool is_constant = true,
              const std::string &func_name = HOUSEPI_NA_STRING);

    /* Set value to variable. */
    void set_value(double x) { (this->*_do_set_value)(x); }

    /* Recover buffered value. */
    void recover_value() { _value = _buffer_value; }

    /* Restore user defined variable initialization. */
    void user_reset()
    {
        _value = _user_value;
        _prop_ptr->user_reset();
    }

    /* Assign a scalar argument to function. */
    void assign_func_attr(const std::string &a, bool is_arg, double x) { _func_ptr->assign_attr(a, is_arg, x); }
    /* Assign a pointer argument to function. */
    void assign_func_attr(const std::string &a, bool is_arg, const double *y) { _func_ptr->assign_attr(a, is_arg, y); }
    /* Assign a function attribute to function. */
    void assign_func_attr(const std::string &a, bool is_arg, std::shared_ptr<Function> f) { _func_ptr->assign_attr(a, is_arg, f); }

    /* Initialize function reshaping. */
    void init_func_reshaping() { _func_ptr->init_reshaping(_is_discrete); }

    //----- Pointer to distribution -----//
    void assign_distr(const std::string &distr_name = HOUSEPI_NA_STRING);

    /* Assign a scalar argument to distribution. */
    void assign_distr_attr(const std::string &a, bool is_arg, double x) { _distr_ptr->assign_attr(a, is_arg, x); }
    /* Assign a pointer argument to distribution. */
    void assign_distr_attr(const std::string &a, bool is_arg, const double *y) { _distr_ptr->assign_attr(a, is_arg, y); }
    /* Assign a pointer argument to distribution. */
    void assign_distr_attr(const std::string &a, bool is_arg, std::shared_ptr<Function> f) { _distr_ptr->assign_attr(a, is_arg, f); }

    void init_distr_reshaping();

    //----- Pointer to proposer -----//
    void assign_proposer(const std::string &rule_name = HOUSEPI_NA_STRING,
                         double rule_scale = HOUSEPI_NA_DOUBLE,
                         double rule_from = HOUSEPI_NA_DOUBLE);

    //===============================//
    //===== METROPOLIS-HASTINGS =====//
    //===============================//

    double get_move_logratio() const { return (this->*_do_move_logratio)(); }
    bool valid_move() const { return (this->*_do_valid_move)(); }

    //----- distr_ptr : prior/model methods -----//
    void rand(rnd_gen &prng) { set_value(_distr_ptr->rand(prng)); }

    double distr_logpdf() const { return _distr_ptr->logpdf(_value); }
    double distr_logratio() const { return _distr_ptr->logratio(_buffer_value, _value); }

    //----- prop_ptr : propose methods -----//
    void propose(rnd_gen &prng) { set_value((this->*_do_propose)(prng)); }
    double propose_logratio() const { return (_prop_ptr->logratio)(_buffer_value, _value); }

    //----- prop_ptr : adaptive methods -----//
    double get_accept_rate() const { return _prop_ptr->get_accept_rate(); }
    void update_accept_rate(bool is_accepted) { _prop_ptr->update_accept_rate(is_accepted); }
    void reset_accept_rate() { _prop_ptr->reset_accept_rate(); }
    void adapt_rule(double ratio) { (_prop_ptr->adapt)(ratio); }

    //================//
    //===== SHOW =====//
    //================//

    const std::string show_name() const { return _variable_name; }
    const std::string show_value() const { return tools::na_to_str(get()); }

    const std::string show_distr_name() const { return _distr_ptr->show_name(); }
    const std::string show_distr_arg_names() const { return _distr_ptr->show_arg_names(); }
    const std::string show_distr_arg_values() const { return _distr_ptr->show_arg_values(); }
    const std::string show_distr_summary() const { return _distr_ptr->show_summary(); }

    const std::string show_min() const { return _distr_ptr->show_min(); }
    const std::string show_max() const { return _distr_ptr->show_max(); }

    const std::string show_func_summary() const { return _func_ptr->show_summary(); }

    const std::string show_prop_name() const { return _prop_ptr->show_name(); }
    const std::string show_prop_arg_names() const { return _prop_ptr->show_arg_names(); }
    const std::string show_prop_arg_values() const { return _prop_ptr->show_arg_values(); }
    const std::string show_prop_scale() const { return _prop_ptr->show_scale(); }

private:
    std::shared_ptr<Function> _func_ptr;      // function
    std::shared_ptr<Distribution> _distr_ptr; // prior\model distribution of the variable
    std::shared_ptr<Proposer> _prop_ptr;      // Monte Carlo proposer

    std::string _variable_name = HOUSEPI_NA_STRING; // name of variable
    double _value = HOUSEPI_NA_DOUBLE;              // value
    double _buffer_value = HOUSEPI_NA_DOUBLE;       // buffered value

    double _user_value = HOUSEPI_NA_DOUBLE; // initial user defined value
    bool _is_discrete = false;              // wheter discrete or not
    bool _is_constant = true;               // wheter constant or not

    //==================//
    //===== GETTER =====//
    //==================//

    double (Variable::*_do_operator_return)(double) const = &Variable::__operator_return_error;
    double __operator_return_error(double ) const { throw std::runtime_error(
        "Variable \"" + show_name() + "\" is not a function."); }
    double __operator_return_func(double x) const { return _func_ptr->return_value(x); }

    double (Variable::*_do_get)() const = &Variable::__get_default;
    double __get_default() const { return _value; }
    double __get_func_return() const { return _func_ptr->return_value(); }

    const double *(Variable::*_do_get_pointer)() const = &Variable::__get_pointer_default;
    const double *__get_pointer_default() const { return &_value; }
    const double *__get_pointer_func_return() const { throw std::runtime_error(
        "Called pointer to function \"" + show_name() + "\"."); }

    //==================//
    //===== SETTER =====//
    //==================//

    void (Variable::*_do_set_value)(double) = &Variable::__set_value_default;
    /* Store the buffer and change value. */
    void __set_value_default(double x)
    {
        _buffer_value = _value;
        _value = x;
    }
    /* Set value to function error. */
    void __set_value_func_error(double) { throw std::runtime_error(
        "Value of function  \"" + show_name() + "\" cannot be changed."); }

    //===============================//
    //===== METROPOLIS-HASTINGS =====//
    //===============================//

    //----- Valid move : functional methods -----//
    bool (Variable::*_do_valid_move)() const = &Variable::__valid_move_default;
    bool __valid_move_default() const { return _distr_ptr->valid(_value); }
    bool __valid_move_draw() const { return true; }

    //----- Move logratio : functional methods -----//
    double (Variable::*_do_move_logratio)() const = &Variable::__move_logratio_default_;
    double __move_logratio_default_() const { return distr_logratio() - propose_logratio(); }
    double __move_logratio_draw() const { return 0.; }

    //----- Propose : functional methods -----//
    double (Variable::*_do_propose)(rnd_gen &) = &Variable::__propose_error;
    double __propose_error(rnd_gen &) { throw std::invalid_argument(
        "Called unknown proposer for variable \"" + show_name() + "\"."); }
    double __propose_draw_sample(rnd_gen &prng) { return (_prop_ptr->propose)(prng, 0); }
    double __propose_from_previous(rnd_gen &prng) { return (_prop_ptr->propose)(prng, _value); }
    double __propose_from_zero(rnd_gen &prng) { return (_prop_ptr->propose)(prng, _prop_ptr->get_from_value()); }
};

#endif // HOUSEPI__VARIABLE__H
