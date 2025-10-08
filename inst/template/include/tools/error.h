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
 * tools/error.h
 *
 */

#ifndef HOUSEPI__TOOLS__ERROR__H
#define HOUSEPI__TOOLS__ERROR__H

#include <string>

namespace hsp
{
    /// Display a single-line warning message.
    void warning(const std::string &txt);
    /// Display a single-line error message.
    void error(const std::string &txt);
}

#endif // HOUSEPI__TOOLS__ERROR__H
