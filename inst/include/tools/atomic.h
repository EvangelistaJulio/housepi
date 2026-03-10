/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <stdexcept>
#include <string>
#include <atomic>

#include "tools/openmp.h"

namespace hsp::atomic
{
    struct Error
    {
    public:
        /// Restore error of the critical scope.
        void restore() noexcept
        {
            has_error_.exchange(false);
            message_.clear();
        }
        /// Thread-safe set for first error encountered.
        void set(const std::string &message)
        {
            if (!has_error_.exchange(true, std::memory_order_relaxed))
            {
#pragma omp critical
                message_ = message;
            }
        }
        /// Overload for exception.
        void set(const std::exception &e) { set(e.what()); }
        /// Wheter it should be skipped or not.
        bool is_pending() const noexcept { return has_error_.load(std::memory_order_relaxed); }
        /// If an error was caught in critical scope, throw an exception.
        void throw_exception()
        {
            if (!has_error_)
                return;
            std::string sms = message_;
            restore();
            throw std::runtime_error(sms);
        }

    private:
        std::atomic<bool> has_error_{false}; /// Wheter error was raised during multi-threading
        std::string message_ = {};           /// Error message for last exception during multi-threading
    };
}
