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
    Proposer(const std::string &distr_name,
             double rule_scale = HOUSEPI_NA_DOUBLE,
             double rule_start_value = HOUSEPI_NA_DOUBLE,
             const std::vector<ushort> &args_to_adapt = {},
             bool is_discrete = false);
    Proposer(std::shared_ptr<Distribution>,
             double rule_scale = HOUSEPI_NA_DOUBLE,
             double rule_start_value = HOUSEPI_NA_DOUBLE,
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
    double __propose_undefined(rnd_gen &prng, double start);
    double __propose_sum_to_previous(rnd_gen &prng, double start);
    double __propose_times_previous(rnd_gen &prng, double start);
    double __propose_draw_rand(rnd_gen &prng, double start);

    //===== logratio : functional methods =====//
    double (Proposer::*_do_logratio)(double, double) const = &Proposer::__logratio_undefined;
    double __logratio_undefined(double xold, double xnew) const;
    double __logratio_pdf(double xold, double xnew) const;
    double __logratio_sum_to_zero(double xold, double xnew) const;
    double __logratio_times_zero(double xold, double xnew) const;
    double __logratio_symm(double xold, double xnew) const;
    double __logratio_lognorm_symm(double xold, double xnew) const;
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
             const std::string &func_name = HOUSEPI_NA_STRING);
    ~Variable();

    //==================//
    //===== GETTER =====//
    //==================//

    double get() const;
    const double *get_pointer() const;
    std::shared_ptr<Function> get_function() const;
    double get_from_value() const;
    bool is_function() const;
    bool is_constant() const;
    bool is_discrete() const;
    bool is_from_zero() const;

    //==================//
    //===== SETTER =====//
    //==================//

    void user_reset();
    void set_value(double x);
    void recover_value();
    void init(const std::string &name, double value, bool is_discrete,
              bool is_constant = true,
              const std::string &func_name = HOUSEPI_NA_STRING);

    void assign_func_attr(const std::string &arg_name, bool is_arg, double x);
    void assign_func_attr(const std::string &arg_name, bool is_arg, const double *y);
    void assign_func_attr(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f);

    //----- Pointer to distribution -----//
    void assign_distr(const std::string &distr_name = HOUSEPI_NA_STRING);
    void assign_distr_attr(const std::string &arg_name, bool is_arg, double x);
    void assign_distr_attr(const std::string &arg_name, bool is_arg, const double *y);
    void assign_distr_attr(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f);
    void init_distr_reshaping();

    //----- Pointer to proposer -----//
    void assign_proposer(const std::string &rule_name = HOUSEPI_NA_STRING,
                         double rule_scale = HOUSEPI_NA_DOUBLE,
                         double rule_from = HOUSEPI_NA_DOUBLE);

    //===============================//
    //===== METROPOLIS-HASTINGS =====//
    //===============================//

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

    //================//
    //===== SHOW =====//
    //================//

    const std::string show_name() const;
    const std::string show_value() const;

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

    double (Variable::*_do_get)() const = &Variable::__get_default;
    double __get_default() const;
    double __get_func_return() const;

    const double *(Variable::*_do_get_pointer)() const = &Variable::__get_pointer_default;
    const double *__get_pointer_default() const;
    const double *__get_pointer_func_return() const;

    //==================//
    //===== SETTER =====//
    //==================//

    void (Variable::*_do_set_value)(double) = &Variable::__set_value_default;
    void __set_value_default(double x);
    void __set_value_func_error(double x);

    //===============================//
    //===== METROPOLIS-HASTINGS =====//
    //===============================//

    //----- Valid move : functional methods -----//
    bool (Variable::*_do_valid_move)() const = &Variable::__valid_move_default;
    bool __valid_move_default() const;
    bool __valid_move_draw() const;

    //----- Move logratio : functional methods -----//
    double (Variable::*_do_move_logratio)() const = &Variable::__move_logratio_default_;
    double __move_logratio_default_() const;
    double __move_logratio_draw() const;

    //----- Propose : functional methods -----//
    double (Variable::*_do_propose)(rnd_gen &) = &Variable::__propose_error;
    double __propose_error(rnd_gen &prng);
    double __propose_draw_sample(rnd_gen &prng);
    double __propose_from_previous(rnd_gen &prng);
    double __propose_from_zero(rnd_gen &prng);
};

#endif // HOUSEPI__VARIABLE__H
