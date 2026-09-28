/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <utility>

#ifndef _STANDALONE
#include <RcppParallel.h>
#include <tbb/task_arena.h>
#endif

namespace housepi::parallel
{
    inline constexpr std::size_t global_workers = 0; /// Worker count deferring to the global budget.

    namespace detail
    {
        /// Run the serial loop, shared by standalone and the single-worker shortcut.
        template <class F>
        inline void serial_for_each(std::size_t n, F &body)
        {
            for (std::size_t i = 0; i < n; ++i)
                body(i);
        }

        /// Run the serial reduction, shared by standalone and the single-worker shortcut.
        /// @note `T` value-initializes to the identity and accumulates via `operator+=`.
        template <class T, class F>
        inline T serial_for_reduce(std::size_t n, F &body)
        {
            T acc{};
            for (std::size_t i = 0; i < n; ++i)
                acc += body(i);
            return acc;
        }
    }

#ifdef _STANDALONE

    inline std::size_t current_thread_index() noexcept { return 0; }

    inline bool in_parallel_region() noexcept { return false; }

    template <class F>
    inline void for_each(std::size_t n, F &&body, std::size_t = global_workers)
    {
        detail::serial_for_each(n, body);
    }

    template <class T, class F>
    inline T for_reduce(std::size_t n, F &&body, std::size_t = global_workers)
    {
        return detail::serial_for_reduce<T>(n, body);
    }

    template <class F>
    inline double for_reduce_sum(std::size_t n, F &&body, std::size_t = global_workers)
    {
        return for_reduce<double>(n, std::forward<F>(body));
    }

#else

    namespace detail
    {
        /// Return number of threads requested.
        /// @note Held in a function-local static.
        inline int &num_threads_ref() noexcept
        {
            static int n = -1;
            return n;
        }

        /// Whether the calling thread is running the body of a parallel dispatch.
        inline bool &region_flag() noexcept
        {
            static thread_local bool inside = false;
            return inside;
        }

        /// Hold the region flag for the lifetime of a worker's range.
        struct RegionGuard
        {
            bool previous;
            RegionGuard() noexcept : previous(region_flag()) { region_flag() = true; }
            ~RegionGuard() { region_flag() = previous; }
        };

        /// Return workers a dispatch may use.
        ///
        /// - own budget: that budget;
        /// - no own budget, inside a parallel region: one worker;
        /// - no own budget, outside: the global budget.
        inline int resolve_workers(std::size_t workers) noexcept
        {
            if (workers != global_workers)
                return static_cast<int>(workers);
            return region_flag() ? 1 : num_threads_ref();
        }

        template <class F>
        struct EachWorker : public RcppParallel::Worker
        {
            F &body;
            explicit EachWorker(F &b) : body(b) {}
            void operator()(std::size_t begin, std::size_t end) override
            {
                RegionGuard guard;
                for (std::size_t i = begin; i < end; ++i)
                    body(i);
            }
        };

        template <class T, class F>
        struct Reducer : public RcppParallel::Worker
        {
            F &body;
            T value{};
            explicit Reducer(F &b) : body(b) {}
            Reducer(const Reducer &other, RcppParallel::Split) noexcept
                : body(other.body) {}
            void operator()(std::size_t begin, std::size_t end) override
            {
                RegionGuard guard;
                T local{};
                for (std::size_t i = begin; i < end; ++i)
                    local += body(i);
                value += local;
            }
            void join(const Reducer &other) noexcept { value += other.value; }
        };
    }

    /// Return maximum number of available threads.
    inline std::size_t max_num_threads() noexcept
    {
        int n = tbb::this_task_arena::max_concurrency();
        return n <= 0 ? 1 : static_cast<std::size_t>(n);
    }

    /// Set number of threads, 0 read as 1.
    inline void set_num_threads(std::size_t n) noexcept
    {
        detail::num_threads_ref() = (n == 0) ? 1 : static_cast<int>(n);
    }

    /// Whether the caller runs inside the body of a parallel dispatch.
    inline bool in_parallel_region() noexcept { return detail::region_flag(); }

    /// Return index of the calling worker, 0 outside any parallel region.
    inline std::size_t current_thread_index() noexcept
    {
        if (detail::num_threads_ref() == 1)
            return 0;
        int idx = tbb::this_task_arena::current_thread_index();
        return idx < 0 ? 0 : static_cast<std::size_t>(idx);
    }

    /// Run `body(i)` for i in [0, n) across `workers` RcppParallel workers.
    /// - single worker: serial loop.
    template <class F>
    inline void for_each(std::size_t n, F &&body, std::size_t workers = global_workers)
    {
        int count = detail::resolve_workers(workers);
        if (count == 1)
            return detail::serial_for_each(n, body);
        detail::EachWorker<F> w(body);
        RcppParallel::parallelFor(0, n, w, 1, count);
    }

    /// Run `body(i)` for i in [0, n) across `workers` workers, reducing into `T` via `operator+=`.
    template <class T, class F>
    inline T for_reduce(std::size_t n, F &&body, std::size_t workers = global_workers)
    {
        int count = detail::resolve_workers(workers);
        if (count == 1)
            return detail::serial_for_reduce<T>(n, body);
        detail::Reducer<T, F> r(body);
        RcppParallel::parallelReduce(0, n, r, 1, count);
        return r.value;
    }

    /// Run `body(i)` for i in [0, n) and return sum of contributions.
    template <class F>
    inline double for_reduce_sum(std::size_t n, F &&body, std::size_t workers = global_workers)
    {
        return for_reduce<double>(n, std::forward<F>(body), workers);
    }
#endif
}
