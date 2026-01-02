/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Author: Francesco Camaglia
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * tools/types.h
 *
 */

#ifndef HOUSEPI__TOOLS__TYPES__H
#define HOUSEPI__TOOLS__TYPES__H

/// DEVELOPER: flag used to develop Rcpp chunks
//#define _HOUSEPI_RCPP

#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>

using int_vec = std::vector<int>;
using dbl_vec = std::vector<double>;
using str_vec = std::vector<std::string>;

template <typename T>
using matrix = std::vector<std::vector<T>>;

// unordered dictionaries
template <typename T>
using u_dict = std::unordered_map<std::string, T>;

#endif // HOUSEPI__TOOLS__TYPES__H
