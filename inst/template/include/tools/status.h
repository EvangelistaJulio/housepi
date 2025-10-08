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
 * tools/status.h
 *
 */

#ifndef HOUSEPI__TOOLS__STATUS__H
#define HOUSEPI__TOOLS__STATUS__H

namespace hsp
{
    namespace status
    {
        constexpr int INFECTED = 1; /// Infected status (>=)
        constexpr int HEALTHY = 0;  /// Healthy status (==)
        constexpr int MISSING = -1; /// Missing status (<=)

        inline bool is_infected(int s) { return s >= INFECTED; }
        inline bool is_healthy(int s) { return s == HEALTHY; }
        inline bool is_missing(int s) { return s <= MISSING; }
    }
} // namespace hsp

#endif // HOUSEPI__TOOLS__STATUS__H
