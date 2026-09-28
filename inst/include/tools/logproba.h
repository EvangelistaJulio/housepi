/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cmath>
#include <stdexcept>

#include "tools/nan.h"

namespace housepi
{
    /// Hold a log-probability and its previous value, for Metropolis moves.
    struct LogProba
    {
        /// Set value to `new_ll`, keeping the previous one.
        /// - NaN `new_ll`: error.
        void set(double new_ll)
        {
            old_ = value_;
            value_ = new_ll;

            if (na::is(new_ll))
                throw std::runtime_error(
                    "Log-probability computation returned NaN value.");
        }
        /// Increment value by `delta_ll`, keeping the previous one.
        /// - NaN result: error.
        void increment(double delta_ll) { set(value_ + delta_ll); }
        /// Restore previous value.
        void recover() noexcept { value_ = old_; }

        /// Return current value.
        double value() const noexcept { return value_; }
        /// Return log-ratio of current to previous value.
        double ratio() const noexcept { return value_ - old_; }

        /// Whether value is infinite.
        bool is_inf() const noexcept { return std::isinf(value_); }

    private:
        double value_ = 0.; /// Current log-probability.
        double old_ = 0.;   /// Previous log-probability.
    };
}
