/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cmath>
#include <stdexcept>

#include "tools/nan.h"

namespace hsp
{
    struct LogLike
    {
        void recover() noexcept { value_ = old_; }
        double log_ratio() const noexcept { return value_ - old_; }
        double value() const noexcept { return value_; }

        /// Set new log-likelihood value. Store old value.
        /// Raise exception if new value is a NaN.
        void set(double new_ll)
        {
            old_ = value_;
            value_ = new_ll;

            if (hsp::na::is(new_ll))
                throw std::runtime_error(
                    "Log-likelihood computation returned NaN value.");
        }
        /// Increment log-likelihood by a delta. Store old value.
        /// Raise exception if new value is a NaN.
        void increment(double delta_ll) { set(value_ + delta_ll); }

        /// Wheter value is infinity.
        bool is_inf() const noexcept { return std::isinf(value_); }

    protected:
        double value_ = 0.; /// Current log-likelihood
        double old_ = 0.;   /// Old log-likelihood
    };
}
