/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__RANDOM__H
#define HOUSEPI__TOOLS__RANDOM__H

#include <random>

namespace hsp
{
#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
    /// Mersenne Twister 19937 pseudo-random number generator for 64bit architecture.
    using rnd_gen = std::mt19937_64;
#else
    /// Mersenne Twister 19937 pseudo-random number generator for 32bit architecture.
    using rnd_gen = std::mt19937;
#endif
}

#endif // HOUSEPI__TOOLS__RANDOM__H