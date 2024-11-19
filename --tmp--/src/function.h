/*
 * function.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef FUNCTION__H
#define FUNCTION__H

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

// Forward declaration of implementation class
class Base_Function_Impl;

// FIXME :
// for the moment it is not possible to provide default pointer argumets
// a solution could be point to internal, at that point the whole pointer thing could be collapsed

/* Parent class to all function classes where the arguements are stored in arrays. */
class _Base_Function_
{
public:
    //=====Constructors=====//
    _Base_Function_();
    virtual ~_Base_Function_();

    //----- FUNCTION NAME -----//
    const std::string &get_name() const;

    //----- ARGUMENTS -----//
    double get_arg(int i) const;
    void set_arg(int i, double new_value);
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;

    //----- RANGE -----//
    const std::string show_range_min() const;
    const std::string show_range_max() const;

    //----- VALID -----//
    bool valid(double x) const;

protected:
    // NOTE : to provide a default order is the only way  I found to enforce an order
    // this is necessary if I want to access arguments by index and not by name (faster...?)
    // FIXME : 1) check if by index faster; 2) collapse overloaded constructors
    // FIXME : for the moment it is not possible to provide default pointer argumets
    // a solution could be point to internal, at that point the whole pointer thing could be collapsed
    void init(const std::string &name,
              const dictargs &args_val,
              const std::vector<std::string> &order,
              const dictargs &domain = {});
    void init(const std::string &name,
              const dictargs_ptr &args_ptr,
              const std::vector<std::string> &order,
              const dictargs &domain = {});

    //----- ARGUMENTS -----//
    // note: just for initial non-default arguments
    void assign_args(const dictargs &new_args_val);
    void assign_args(const dictargs_ptr &new_args_ptr);

    //----- RANGE -----//
    ushort set_range(const dictargs &new_range);
    double get_range_min() const;
    double get_range_max() const;

private:
    //===== Pointer to Implementation =====//
    std::unique_ptr<Base_Function_Impl> pImpl;
};

#endif // FUNCTION__H