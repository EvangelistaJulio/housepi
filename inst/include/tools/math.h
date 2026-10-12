/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cmath>

//////  Constants  //////

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288 /// Pi.
#endif

#ifndef M_LN2
#define M_LN2 0.693147180559945309417232121458 /// log(2).
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.41421356237309504880168872420969808 /// sqrt(2).
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362105 /// 1 / sqrt(2).
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI 0.398942280401432677939946059934 /// 1 / sqrt(2 pi).
#endif

#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 /// log(sqrt(2 pi)).
#endif

namespace housepi::math
{
    /// Return `exp(x) - 1`, to a relative precision of 2^-42 or better.
    ///
    /// - `|x| < 2^-10`: `std::expm1`, where `exp(x) - 1` cancels;
    /// - otherwise: `exp(x) - 1`, cheaper than `std::expm1`.
    inline double expm1(double x) { return std::fabs(x) < 0x1p-10 ? std::expm1(x) : std::exp(x) - 1.; }

    /// Return `log(1 - exp(-r))` for `r >= 0`, to a relative precision of 2^-42 or better.
    ///
    /// Split at `log(2)` after @cite machler_accurately_2012:
    /// - `r < log(2)`: `log(-expm1(-r))`, where `1 - exp(-r)` cancels;
    /// - otherwise: `log1p(-exp(-r))`, where `log` of a value near 1 cancels.
    inline double log1mexp(double r) { return r < M_LN2 ? std::log(-expm1(-r)) : std::log1p(-std::exp(-r)); }
}
