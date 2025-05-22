/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * function.h
 *
 */

#ifndef HOUSEPI__FUNCTION__H
#define HOUSEPI__FUNCTION__H

#include <string>
#include <vector>

#include "tools.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

// forward declaration
class Function;

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Attribute : class ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Attribute parent class. */
class Attribute
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Attribute() {};
    virtual ~Attribute() = default;

    //==================//
    //===== GETTER =====//
    //==================//

    double get() const;
    const double *get_pointer() const;
    const std::shared_ptr<Function> get_function() const;
    bool is_scalar() const;
    bool is_pointer() const;
    bool is_function() const;
    const std::string show() const;
    const std::string show_name() const;

    //==================//
    //===== SETTER =====//
    //==================//

    void set(double x);
    void set_name(const std::string &name);
    void assign_scalar(double x);
    void assign_pointer(const double *y);
    void assign_function(const std::shared_ptr<Function> f);

private:
    std::string _attribute_name = HOUSEPI_NA_STRING; // name of argument
    double _value = HOUSEPI_NA_DOUBLE;              // internal value
    const double *_ptr_to_value = nullptr;          // pointer to value
    std::shared_ptr<Function> _func_ptr = nullptr;  // pointer to function

    // Functional getter //
    double (Attribute::*_do_get)() const = &Attribute::__get_default;
    double __get_default() const;
    double __get_from_pointer() const;
    double __get_from_function() const;

    // Functional pointer getter //
    const double * (Attribute::*_do_get_pointer)() const = &Attribute::__get_pointer_default;
    const double * __get_pointer_default() const;
    const double * __get_pointer_from_pointer() const;
    const double * __get_pointer_error() const;

    // Functional setter //
    void (Attribute::*_do_set)(double) = &Attribute::__set_default;
    void __set_default(double x);
    void __set_error(double x);
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : class //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Function parent class. Arguements stored in arrays. */
class Function
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Function() {};
    virtual ~Function() = default;

    //==================//
    //===== GETTER =====//
    //==================//

    virtual double return_value(const dbl_vec &input) const;

    ushort find(const std::string &arg_name) const;
    ushort assert_arg(const std::string &arg_name) const;
    double get_arg(ushort i) const;
    double get_min() const;
    double get_max() const;
    const double *get_ptr_to_arg(ushort i) const;
    const double *get_ptr_to_min() const;
    const double *get_ptr_to_max() const;
    bool below(double x) const;
    bool above(double x) const;
    bool valid(double x) const;
    ushort which_valid_case() const;

    //=============================//
    //===== ARGUMENTS & RANGE =====//
    //=============================//
    void assign_attr_scalar(const std::string &arg_name, bool is_arg, double value);
    void assign_attr_pointer(const std::string &arg_name, bool is_arg, const double *ptr_to_value);
    void assign_attr_function(const std::string &arg_name, bool is_arg, std::shared_ptr<Function> f);
    void set_arg(ushort i, double new_value);
    void init_uniform_range();
    void init_functional_range(bool is_min_open, bool is_max_open);

    //===== SHOW =====//
    const std::string show_name() const;
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_min() const;
    const std::string show_max() const;
    const std::string show_range() const;
    const std::string show_summary() const;

protected:
    void init_defaults(const std::string &name,
                       const str_vec &arg_name,
                       const dbl_vec &arg_value,
                       double domain_min,
                       double domain_max,
                       bool is_min_open,
                       bool is_max_open);

private:
    std::string _function_name = HOUSEPI_NA_STRING;

    //===== Arguments =====//
    Attribute *_arguments = nullptr;
    ushort _n_args = 0;

    //===== Domain =====//
    Attribute _range_min = Attribute();
    Attribute _range_max = Attribute();

    //=============================//
    //===== ARGUMENTS & RANGE =====//
    //=============================//
    void _assign_arg_scalar(const std::string &arg_name, double value);
    void _assign_arg_pointer(const std::string &arg_name, const double *ptr_to_value);
    void _assign_arg_function(const std::string &arg_name, const std::shared_ptr<Function> f);

    void _assign_range_scalar(const std::string &range_name, double lower);
    void _assign_range_pointer(const std::string &range_name, const double *ptr_to_lower);
    void _assign_range_function(const std::string &range_name, const std::shared_ptr<Function> f);

    //===== VALID =====//
    // note: for the moment range extremes can't be one open and other close

    //---- alternative valid range methods -----//
    bool (Function::*_do_below)(double) const = &Function::__below_never;
    bool (Function::*_do_above)(double) const = &Function::__above_never;
    bool (Function::*_do_valid)(double) const = &Function::__valid_always;

    //----- range : none ------//
    bool __below_never(double x) const
    {
        UNUSED(x);
        return false;
    }
    bool __above_never(double x) const
    {
        UNUSED(x);
        return false;
    }
    bool __valid_always(double x) const
    {
        UNUSED(x);
        return true;
    }

    //----- range : minimum -----//
    bool __below_open_min(double x) const { return (x <= _range_min.get()); }
    bool __below_close_min(double x) const { return (x < _range_min.get()); }
    bool __valid_open_min(double x) const { return (x > _range_min.get()); }
    bool __valid_close_min(double x) const { return (x >= _range_min.get()); }

    //----- range : maximum -----//
    bool __above_open_max(double x) const { return (x >= _range_max.get()); }
    bool __above_close_max(double x) const { return (x > _range_max.get()); }
    bool __valid_open_max(double x) const { return (x < _range_max.get()); }
    bool __valid_close_max(double x) const { return (x <= _range_max.get()); }

    //----- range : both -----//
    bool __valid_both_open(double x) const { return (x > _range_min.get()) && (x < _range_max.get()); }
    bool __valid_open_min_close_max(double x) const { return (x > _range_min.get()) && (x <= _range_max.get()); }
    bool __valid_close_min_open_max(double x) const { return (x >= _range_min.get()) && (x < _range_max.get()); }
    bool __valid_both_close(double x) const { return (x >= _range_min.get()) && (x <= _range_max.get()); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  function switch  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Function with default scalar arguments. */
std::shared_ptr<Function> function_switch(const std::string &name = HOUSEPI_NA_STRING);

#endif // HOUSEPI__FUNCTION__H
