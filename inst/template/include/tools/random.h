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
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * tools/random.h
 *
 */

#ifndef HOUSEPI__TOOLS__RANDOM__H
#define HOUSEPI__TOOLS__RANDOM__H

#include <random>

// Use pseudo-random number generator for 64bit architecture.
#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
using rnd_gen = std::mt19937_64;
#else
using rnd_gen = std::mt19937;
#endif

#endif // HOUSEPI__TOOLS__RANDOM__H