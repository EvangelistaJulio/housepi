/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <utility>

#ifndef _STANDALONE
#include <RcppParallel.h>
#endif

namespace housepi::parallel
{
    namespace detail
    {
        /// Serial fallback loop, shared by standalone and the is_serial() shortcut.
        template <class F>
        inline void serial_for_each(std::size_t n, F &body)
        {
            for (std::size_t i = 0; i < n; ++i)
                body(i);
        }

        /// Serial sum-reduction, shared by standalone and the is_serial() shortcut.
        template <class F>
        inline double serial_for_reduce_sum(std::size_t n, F &body)
        {
            double acc = 0.0;
            for (std::size_t i = 0; i < n; ++i)
                acc += body(i);
            return acc;
        }
    }

#ifdef _STANDALONE

    inline std::size_t current_thread_index() noexcept { return 0; }

    template <class F>
    inline void for_each(std::size_t n, F &&body) { detail::serial_for_each(n, body); }

    template <class F>
    inline double for_reduce_sum(std::size_t n, F &&body)
    {
        return detail::serial_for_reduce_sum(n, body);
    }

#else

    namespace detail
    {
        /// Number of threads requested.
        /// NOTE: stored at file scope.
        inline int &num_threads_ref() noexcept
        {
            static int n = -1;
            return n;
        }

        /// Whether running on a single thread.
        inline bool is_serial() noexcept { return num_threads_ref() == 1; }

        template <class F>
        struct EachWorker : public RcppParallel::Worker
        {
            F &body;
            explicit EachWorker(F &b) : body(b) {}
            void operator()(std::size_t begin, std::size_t end) override
            {
                for (std::size_t i = begin; i < end; ++i)
                    body(i);
            }
        };

        template <class F>
        struct SumReducer : public RcppParallel::Worker
        {
            F &body;
            double value = 0.0;
            explicit SumReducer(F &b) : body(b) {}
            SumReducer(const SumReducer &other, RcppParallel::Split) noexcept
                : body(other.body) {}
            void operator()(std::size_t begin, std::size_t end) override
            {
                double local = 0.0;
                for (std::size_t i = begin; i < end; ++i)
                    local += body(i);
                value += local;
            }
            void join(const SumReducer &other) noexcept { value += other.value; }
        };
    }

    /// Return maximum number of available threads.
    inline std::size_t max_num_threads() noexcept
    {
        int n = tbb::this_task_arena::max_concurrency();
        return n <= 0 ? 1 : static_cast<std::size_t>(n);
    }

    /// Set number of threads.
    inline void set_num_threads(std::size_t n) noexcept
    {
        detail::num_threads_ref() = (n == 0) ? 1 : static_cast<int>(n);
    }

    /// Index of the calling worker.
    /// Returns 0 when called outside any parallel region.
    inline std::size_t current_thread_index() noexcept
    {
        if (detail::is_serial())
            return 0;
        int idx = tbb::this_task_arena::current_thread_index();
        return idx < 0 ? 0 : static_cast<std::size_t>(idx);
    }

    /// Run `body(i)` for i in [0, n) across RcppParallel workers.
    template <class F>
    inline void for_each(std::size_t n, F &&body)
    {
        if (detail::is_serial())
            return detail::serial_for_each(n, body);
        detail::EachWorker<F> w(body);
        RcppParallel::parallelFor(0, n, w, 1, detail::num_threads_ref());
    }

    /// Run `body(i)` for i in [0, n) and return sum of contributions.
    template <class F>
    inline double for_reduce_sum(std::size_t n, F &&body)
    {
        if (detail::is_serial())
            return detail::serial_for_reduce_sum(n, body);
        detail::SumReducer<F> r(body);
        RcppParallel::parallelReduce(0, n, r, 1, detail::num_threads_ref());
        return r.value;
    }
#endif
}
