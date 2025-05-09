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

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Argument : class ////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Argument
{
public:
    //=======================//
    //===== CONSTRUCTOR =====//
    //=======================//

    Argument() {};
    virtual ~Argument() = default;

    //==================//
    //===== GETTER =====//
    //==================//

    double get() const;
    const double *get_ptr() const;
    bool is_internal_ptr() const;
    const std::string show() const;
    const std::string show_name() const;

    //==================//
    //===== SETTER =====//
    //==================//

    void set(double x);
    void set_name(const std::string &name);
    void assign_ptr(const double *y);
    void restore_internal_ptr();

private:
    std::string _argument_name = HHEPI_NA_STRING;
    double _value = HHEPI_NA_DOUBLE;
    const double *_ptr_to_value = &_value;

    void (Argument::*_func_set_arg)(double) = &Argument::__set_arg_default;
    void __set_arg_default(double x);
    void __set_arg_pointer(double x);
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : class //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Parent class to all function classes where the arguements are stored in arrays. */
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

    //===== ARGUMENTS =====//
    void assign_scalar_arg(const std::string &arg_name, double value);
    void assign_pointer_arg(const std::string &arg_name, const double *ptr_to_value);
    void set_arg(ushort i, double new_value);

    //===== RANGE =====//
    void assign_scalar_range_min(double lower);
    void assign_pointer_range_min(const double *ptr_to_lower);
    void assign_scalar_range_max(double upper);
    void assign_pointer_range_max(const double *ptr_to_upper);
    void init_functional_range(bool is_min_open, bool is_max_open);

    //===== SHOW =====//
    const std::string show_name() const;
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_min() const;
    const std::string show_max() const;
    const std::string show_range() const;

protected:
    void init_defaults(const std::string &name,
                       const str_vec &arg_name,
                       const dbl_vec &arg_value,
                       double domain_min,
                       double domain_max,
                       bool is_min_open,
                       bool is_max_open);

private:
    std::string _function_name = HHEPI_NA_STRING;

    //===== Arguments =====//
    Argument *_arguments = nullptr;
    ushort _n_args = 0;

    //===== Domain =====//
    Argument _range_min = Argument();
    Argument _range_max = Argument();

    //===== VALID =====//
    // note: for the moment range extremes can't be one open and other close

    //---- alternative valid range methods -----//
    bool (Function::*_func_below)(double) const = &Function::__below_never;
    bool (Function::*_func_above)(double) const = &Function::__above_never;
    bool (Function::*_func_valid)(double) const = &Function::__valid_always;

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

#endif // HOUSEPI__FUNCTION__H
