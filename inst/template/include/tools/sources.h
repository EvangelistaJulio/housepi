/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * tools/sources.h
 *
 */

#ifndef HOUSEPI__TOOLS__SOURCES__H
#define HOUSEPI__TOOLS__SOURCES__H

#include <string>

namespace hsp
{
    // Source checking functions.
    namespace sources
    {
        inline bool is_augmented(const std::string &s) { return s == "A"; }
        inline bool is_data(const std::string &s) { return s == "D"; }
        inline bool is_group_rate(const std::string &s) { return s == "G"; }
        inline bool is_numeric(const std::string &s) { return s == "N"; }
        inline bool is_parameter(const std::string &s) { return s == "P"; }
    }
}

#endif // HOUSEPI__TOOLS__SOURCES__H
