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

    template <class F>
    inline void for_each(std::size_t n, F &&body, std::size_t /*workers*/)
    {
        detail::serial_for_each(n, body);
    }

    template <class T, class F>
    inline T for_reduce(std::size_t n, F &&body, std::size_t /*workers*/)
    {
        return detail::serial_for_reduce<T>(n, body);
    }

    template <class F>
    inline double for_reduce_sum(std::size_t n, F &&body, std::size_t workers)
    {
        return for_reduce<double>(n, std::forward<F>(body), workers);
    }

#else

    namespace detail
    {
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

    /// Return index of the calling worker in its task arena, 0 outside any arena.
    /// @note Each dispatch runs in its own arena of `workers` slots: the index stays below `workers`.
    inline std::size_t current_thread_index() noexcept
    {
        int idx = tbb::this_task_arena::current_thread_index();
        return idx < 0 ? 0 : static_cast<std::size_t>(idx);
    }

    /// Run `body(i)` for i in [0, n) across `workers` RcppParallel workers.
    /// - at most one worker: serial loop on the calling thread.
    template <class F>
    inline void for_each(std::size_t n, F &&body, std::size_t workers)
    {
        if (workers <= 1)
            return detail::serial_for_each(n, body);
        detail::EachWorker<F> w(body);
        RcppParallel::parallelFor(0, n, w, 1, static_cast<int>(workers));
    }

    /// Run `body(i)` for i in [0, n) across `workers` workers, reducing into `T` via `operator+=`.
    /// - at most one worker: serial loop on the calling thread.
    template <class T, class F>
    inline T for_reduce(std::size_t n, F &&body, std::size_t workers)
    {
        if (workers <= 1)
            return detail::serial_for_reduce<T>(n, body);
        detail::Reducer<T, F> r(body);
        RcppParallel::parallelReduce(0, n, r, 1, static_cast<int>(workers));
        return r.value;
    }

    /// Run `body(i)` for i in [0, n) and return sum of contributions.
    template <class F>
    inline double for_reduce_sum(std::size_t n, F &&body, std::size_t workers)
    {
        return for_reduce<double>(n, std::forward<F>(body), workers);
    }
#endif
}
