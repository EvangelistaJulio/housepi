/*
 * function.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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
    std::string name = NULL_STRING;

    //==================//
    //===== GETTER =====//
    //==================//

    double get() const;
    const double *get_ptr() const;
    bool is_internal_ptr() const;
    const std::string show() const;

    //==================//
    //===== SETTER =====//
    //==================//

    void set(double x);
    void assign_ptr(const double *y);
    void restore_internal_ptr();

private:
    double value = NULL_DOUBLE;
    const double *ptr_to_value = &value;

    void (Argument::*_func_set_)(double) = &Argument::_set_arg_default_;
    void _set_arg_default_(double x);
    void _set_arg_pointer_(double x);
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Base Function : class //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Parent class to all function classes where the arguements are stored in arrays. */
class Function
{
public:
    //===== CONSTRUCTOR =====//
    Function();
    virtual ~Function();

    //===== GETTER =====//
    ushort find(const std::string &arg_name) const;
    double get_arg(ushort i) const;
    double get_min() const;
    double get_max() const;
    const double *get_ptr_to_arg(ushort i) const;
    const double *get_ptr_to_min() const;
    const double *get_ptr_to_max() const;
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
    std::string name = NULL_STRING;

    //===== Arguments =====//
    Argument *arguments = nullptr;
    ushort n_args = 0;

    //===== Domain =====//
    Argument range_min = Argument();
    Argument range_max = Argument();

    //===== VALID =====//
    // note: for the moment range extremes can't be one open and other close

    //---- alternative valid range methods -----//
    bool (Function::*_func_valid_range_)(double) const = &Function::_valid_range_always_;
    bool _valid_range_always_(double x) const
    {
        UNUSED(x);
        return true;
    }

    //----- range : minimum -----//
    bool _valid_range_open_min_(double x) const { return (x > range_min.get()); }
    bool _valid_range_close_min_(double x) const { return (x >= range_min.get()); }

    //----- range : maximum -----//
    bool _valid_range_open_max_(double x) const { return (x < range_max.get()); }
    bool _valid_range_close_max_(double x) const { return (x <= range_max.get()); }

    //----- range : both -----//
    bool _valid_range_both_open_(double x) const { return (x > range_min.get()) && (x < range_max.get()); }
    bool _valid_range_open_min_close_max_(double x) const { return (x > range_min.get()) && (x <= range_max.get()); }
    bool _valid_range_close_min_open_max_(double x) const { return (x >= range_min.get()) && (x < range_max.get()); }
    bool _valid_range_both_close_(double x) const { return (x >= range_min.get()) && (x <= range_max.get()); }
};

#endif // HOUSEPI__FUNCTION__H
