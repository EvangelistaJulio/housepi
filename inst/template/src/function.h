/*
 * function.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
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
    double get_arg(int i) const;
    bool valid(double x) const;

    //===== ARGUMENTS =====//
    void assign_args(const dictargs &new_args_val);
    void assign_args(const dictargs_ptr &new_args_ptr);
    void set_arg(int i, double new_value);

    //===== SHOW =====//
    const std::string show_arg_names() const;
    const std::string show_arg_values() const;
    const std::string show_range_min() const;
    const std::string show_range_max() const;

protected:
    // NOTE : to provide a default order is the only way I found to enforce an order
    // this is necessary if I want to access arguments by index and not by name (faster...?)
    // UPDATE : 1) check if by index faster; 2) collapse overloaded constructors
    // UPDATE : for the moment it is not possible to provide default pointer argumets
    // a solution could be point to internal, at that point the whole pointer thing could be collapsed
    void init(const std::string &name,
              const dictargs &args_val,
              const std::vector<std::string> &order,
              const dictargs &domain = {});
    void init(const std::string &name,
              const dictargs_ptr &args_ptr,
              const std::vector<std::string> &order,
              const dictargs &domain = {});

    //----- RANGE -----//
    ushort set_range(const dictargs &new_range);
    double get_range_min() const;
    double get_range_max() const;

private:
    //===== Pointer to Implementation =====//
    std::unique_ptr<Base_Function_Impl> pImpl;
};

#endif // HOUSEPI__FUNCTION__H
