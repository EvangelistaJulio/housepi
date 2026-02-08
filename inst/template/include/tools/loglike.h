/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <limits>
#include <stdexcept>

#include "tools/nan.h"

namespace hsp
{
    struct LogLike
    {
        void recover() noexcept { value_ = buffer_; }
        double logratio() const noexcept { return value_ - buffer_; }
        double value() const noexcept { return value_; }

        /// Set new log-likelihood value. Buffer old value.
        /// Raise exception if new value is a NaN.
        void set(double new_ll)
        {
            if (hsp::na::is(new_ll))
                throw std::runtime_error(
                    "Log-likelihood computation returned NaN value.");

            buffer_ = value_;
            value_ = new_ll;
        }
        /// Increment log-likelihood by a delta. Buffer old value.
        void increment(double delta_ll) { set(value_ + delta_ll); }

        /// Wheter value is infinity.
        bool is_inf() const noexcept { return  std::isinf(value_); }

    protected:
        double value_ = 0.;  /// Current log-likelihood
        double buffer_ = 0.; /// Buffered log-likelihood
    };
}