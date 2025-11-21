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
 * tools/math_combin.h
 *
 */

#ifndef HOUSEPI__TOOLS_MATH_COMBINATORIAL__H
#define HOUSEPI__TOOLS_MATH_COMBINATORIAL__H

#include <cmath>

namespace numr
{
    /// Standard definition for the factorial.
    inline unsigned long factorial(unsigned long n)
    {
        unsigned long result = 1;
        for (unsigned long i = 2; i <= n; ++i)
            result *= i;
        return result;
    }

    /// Binomial coefficient of N elements in k positions
    inline unsigned long binom_coeff(unsigned N, unsigned k)
    {
        if (N == 0 || k == N)
            return 1;
        // note: use symmetry
        if (k > N - k)
            k = N - k;
        unsigned long result = 1;
        for (unsigned i = 1; i <= k; ++i)
            result = result * (N - k + i) / i;
        return result;
    }
}
#endif // HOUSEPI__TOOLS_MATH_COMBINATORIAL__H
