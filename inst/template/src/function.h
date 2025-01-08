/*
 * function.h
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#ifndef HOUSEPI__FUNCTION__H
#define HOUSEPI__FUNCTION__H

#include <memory>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <chrono>

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

/* Forward declaration of implementationl. */
class Base_Function_Impl;

/* Parent class to all function classes where the arguements are stored in arrays. */
class Base_Function
{
public:
    //===== CONSTRUCTOR =====//
    Base_Function();
    virtual ~Base_Function();

    //===== GETTER =====//
    const std::string &get_name() const;
    bool valid(double x) const;
    double get_arg(ushort i) const;
    double get_range_min() const;
    double get_range_max() const;

    //===== ARGUMENTS =====//
    void assign_scalar_arg(const std::string &arg_name, double value);
    void assign_pointer_arg(const std::string &arg_name, const double *ptr_to_value);
    void set_arg(ushort i, double new_value);
    ushort set_range(double lower_value = UNDEFINED_DOUBLE,
                     double upper_value = UNDEFINED_DOUBLE);

    //===== SHOW =====//
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_range_min() const;
    const std::string show_range_max() const;

protected:
    //===== Pointer to Implementation =====//
    std::unique_ptr<Base_Function_Impl> pImpl;

    void _init_defaults_(const std::string &name,
                         const str_vec &arg_name,
                         const dbl_vec &arg_value,
                         double domain_min,
                         double domain_max);
};

#endif // HOUSEPI__FUNCTION__H
