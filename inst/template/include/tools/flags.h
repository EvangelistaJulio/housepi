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
 * tools/flags.h
 *
 */

#ifndef HOUSEPI__TOOLS__FLAGS__H
#define HOUSEPI__TOOLS__FLAGS__H

#include <string>

/// Compile-time list of standard housepi flags.
class CompiledFlagsHousepi
{
public:
    CompiledFlagsHousepi();

    /// Get compilation flags.
    std::string get_flags() { return flags_; }

private:
    std::string flags_;
};

#endif // HOUSEPI__TOOLS__FLAGS__H
