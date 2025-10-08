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
 * tools/params.h
 *
 */

#ifndef HOUSEPI__TOOLS__PARAMS__H
#define HOUSEPI__TOOLS__PARAMS__H

namespace hsp
{
    // Parameter names.
    namespace params
    {
        constexpr const char *X = "x";          /// Evaluation parameter.
        constexpr const char *LOCATION = "loc"; /// Location parameter.
        constexpr const char *MIN = "min";      /// Minimum value.
        constexpr const char *MAX = "max";      /// Maximum value.
    }
}

#endif // HOUSEPI__TOOLS__PARAMS__H
