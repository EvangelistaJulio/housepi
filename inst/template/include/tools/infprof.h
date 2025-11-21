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
 * tools/infprof.h
 *
 */

#ifndef HOUSEPI__TOOLS__INFPROF__H
#define HOUSEPI__TOOLS__INFPROF__H

namespace hsp
{
    namespace infprof
    {
        /// Augmented instantaneous infectivity profile.
        constexpr const char *INSTANTANEOUS = ".profile_inst";
        /// Augmented cumulative infectivity profile.
        constexpr const char *CUMULATIVE = ".profile_cumul";
        /// Default size for lookup table for discrete case
        constexpr unsigned LOOKUP_SIZE = 30;
    }
}

#endif // HOUSEPI__TOOLS__INFPROF__H
