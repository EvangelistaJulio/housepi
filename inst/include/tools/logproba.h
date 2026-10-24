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
        /// * NaN or `+Inf` `new_ll` -> error.
        void set(double new_ll)
        {
            old_ = value_;
            value_ = new_ll;

            // NaN compares false: caught with +Inf
            if (!(new_ll < INFINITY))
                throw_invalid_(new_ll);
        }
        /// Increment value by `delta_ll`, keeping the previous one.
        /// * NaN or `+Inf` result -> error.
        void increment(double delta_ll) { set(value_ + delta_ll); }
        /// Restore previous value.
        void recover() noexcept { value_ = old_; }

        /// Return current value.
        double value() const noexcept { return value_; }
        /// Return log-ratio of current to previous value.
        double ratio() const noexcept { return value_ - old_; }

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

    /// Log-likelihood terms a move recomputes.
    ///
    /// `transm` marks a move reaching the transmission cache, with a meaning per move:
    /// * parameter -> household caches and infectivity profile refreshed;
    /// * latent -> member cache refreshed, transmission recomputed.
    /// @note A parameter move recomputes transmission whatever `transm`.
    struct LogLikTerms
    {
        bool obs = false;    /// Observation.
        bool transm = false; /// Transmission cache, see above.
        bool latent = false; /// Latent model.

        /// Add the terms of `other`.
        LogLikTerms &operator|=(const LogLikTerms &other) noexcept
        {
            obs = obs || other.obs;
            transm = transm || other.transm;
            latent = latent || other.latent;
            return *this;
        }
    };

    /// Observation, transmission and latent-model log-likelihoods, or their deltas.
    struct LogLik
    {
        double obs = 0.;    /// Observation.
        double transm = 0.; /// Transmission.
        double latent = 0.; /// Latent model.

        LogLik &operator+=(const LogLik &other) noexcept
        {
            obs += other.obs;
            transm += other.transm;
            latent += other.latent;
            return *this;
        }
        /// Return every term multiplied by `w`.
        LogLik operator*(double w) const noexcept { return {w * obs, w * transm, w * latent}; }
        /// Return the sum of the terms.
        double value() const noexcept { return obs + transm + latent; }
    };

    /// Log-likelihood totals over households, each with its previous value.
    /// * transmission -> always set, recovered and in the ratio;
    /// * observation, latent model -> when in `terms`.
    struct LogLikTotal
    {
        LogProba obs;    /// Observation.
        LogProba transm; /// Transmission.
        LogProba latent; /// Latent model.

        /// Set the totals to `ll`.
        void set(const LogLik &ll, const LogLikTerms &terms)
        {
            transm.set(ll.transm);
            if (terms.obs)
                obs.set(ll.obs);
            if (terms.latent)
                latent.set(ll.latent);
        }
        /// Increment every total by `delta`.
        void increment(const LogLik &delta)
        {
            obs.increment(delta.obs);
            transm.increment(delta.transm);
            latent.increment(delta.latent);
        }
        /// Restore the previous totals.
        void recover(const LogLikTerms &terms) noexcept
        {
            transm.recover();
            if (terms.obs)
                obs.recover();
            if (terms.latent)
                latent.recover();
        }

        /// Return the sum of the current totals.
        double value() const noexcept { return obs.value() + transm.value() + latent.value(); }
        /// Return the log-ratio of current to previous totals.
        double ratio(const LogLikTerms &terms) const noexcept
        {
            double r = transm.ratio();
            if (terms.obs)
                r += obs.ratio();
            if (terms.latent)
                r += latent.ratio();
            return r;
        }
    };
}
