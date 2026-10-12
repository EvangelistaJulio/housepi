/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cmath>

#include "tools/messages.h"

namespace housepi
{
    /// Hold a log-probability and its previous value, for Metropolis moves.
    struct LogProba
    {
        /// Set value to `new_ll`, keeping the previous one.
        /// - NaN or `+Inf` `new_ll`: error.
        void set(double new_ll)
        {
            old_ = value_;
            value_ = new_ll;

            // NaN compares false: caught with +Inf
            if (!(new_ll < INFINITY))
                throw_invalid_(new_ll);
        }
        /// Increment value by `delta_ll`, keeping the previous one.
        /// - NaN or `+Inf` result: error.
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

        /// Throw the error for an invalid log-probability `ll`, NaN or `+Inf`.
        [[noreturn]] static void throw_invalid_(double ll)
        {
            if (std::isnan(ll))
                msg::stop("Log-likelihood is not a number.");
            msg::stop("Log-likelihood is infinite.");
        }
    };
}
