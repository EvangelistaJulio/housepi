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
 * tools/discrete.h
 *
 */

#ifndef HOUSEPI__TOOLS__DISCR__H
#define HOUSEPI__TOOLS__DISCR__H

namespace hsp
{
    // Discretization types.
    namespace discr
    {
        constexpr unsigned short NONE = 0;    /// None
        constexpr unsigned short LEFT = 1;    /// Left
        constexpr unsigned short CENTER = 2;  /// Center
        constexpr unsigned short DEFINED = 3; /// Defined
    }
}

#endif // HOUSEPI__TOOLS__DISCR__H
