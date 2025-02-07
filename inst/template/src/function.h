/*
 * function.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__FUNCTION__H
#define HOUSEPI__FUNCTION__H

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
    void init_range(bool is_min_open, bool is_max_open);

    //===== SHOW =====//
    const std::string show_name() const;
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_min() const;
    const std::string show_max() const;

protected:
    void init_defaults(const std::string &name,
                         const str_vec &arg_name,
                         const dbl_vec &arg_value,
                         double domain_min,
                         double domain_max,
                         bool is_min_open,
                         bool is_max_open);
private:
    //===== Pointer to Implementation =====//
    class Impl; // forward declaration
    std::unique_ptr<Impl> pImpl;
};

#endif // HOUSEPI__FUNCTION__H
