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
    /// Log-probability used for Metropolis. Aware of old value.
    struct LogProba
    {
        /// Set new log-probability value. Store old value.
        /// Raise exception if new value is a NaN.
        void set(double new_ll)
        {
            old_ = value_;
            value_ = new_ll;

            if (na::is(new_ll))
                throw std::runtime_error(
                    "Log-probability computation returned NaN value.");
        }
        /// Increment log-probability by a delta. Store old value.
        /// Raise exception if new value is a NaN.
        void increment(double delta_ll) { set(value_ + delta_ll); }
        /// Recover old value.
        void recover() noexcept { value_ = old_; }

        /// Return current value of log-probability.
        double value() const noexcept { return value_; }
        /// Return log-ratio current vs old values.
        double ratio() const noexcept { return value_ - old_; }

        /// Wheter value is infinity.
        bool is_inf() const noexcept { return std::isinf(value_); }

    private:
        double value_ = 0.; /// Current log-likelihood
        double old_ = 0.;   /// Old log-likelihood
    };

    /// Accumulator for log-densities computations.
    struct LogLikAccumulator
    {
        double obs_ll = 0.;    /// observation
        double transm_ll = 0.; /// transmission
        double latent_ll = 0.;   /// latent

        LogLikAccumulator &operator+=(const LogLikAccumulator &other) noexcept
        {
            obs_ll += other.obs_ll;
            transm_ll += other.transm_ll;
            latent_ll += other.latent_ll;
            return *this;
        }
    };
}
