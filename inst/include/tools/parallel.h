/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <atomic>
#include <chrono>
#include <cstddef>
#include <exception>
#include <future>
#include <utility>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/task_arena.h>

namespace housepi::parallel
{
    /// First error of a run of units.
    struct Failure
    {
        std::exception_ptr error; /// Exception that stopped the run, null where every unit completed.
        std::size_t index = 0;    /// Unit that raised it.
    };

    /// Thrown by a unit that stops because another unit failed.
    struct Cancelled
    {
    };

    namespace detail
    {
        using Range = tbb::blocked_range<std::size_t>;

        /// Run `f()` on the calling thread, in slot 0 of a new arena of `n_slots` slots.
        /// * at most one slot -> `f()` in the current arena.
        /// @return Value of `f()`.
        template <class F>
        inline auto in_arena(std::size_t n_slots, F &&f)
        {
            // single slot: no arena needed
            if (n_slots <= 1)
                return f();
            return tbb::task_arena(static_cast<int>(n_slots)).execute(f);
        }

        /// Run `f()` in an arena of slots.
        /// * current arena has `n_slots` slots -> `f()` in place;
        /// * otherwise -> `f()` in a new arena.
        /// @return Value of `f()`.
        template <class F>
        inline auto in_arena_of(std::size_t n_slots, F &&f)
        {
            // current arena has `slots`: run in place
            if (tbb::this_task_arena::max_concurrency() == static_cast<int>(n_slots))
                return f();
            return in_arena(n_slots, f);
        }
    }

    /// Return maximum number of available threads.
    inline std::size_t max_num_threads() noexcept
    {
        int n = tbb::this_task_arena::max_concurrency();
        return n <= 0 ? 1 : static_cast<std::size_t>(n);
    }

    /// Return index of the calling worker in its task arena, 0 outside any arena.
    /// @note Each dispatch runs in an arena of `workers` slots: the index stays below `workers`.
    inline std::size_t current_thread_index() noexcept
    {
        int idx = tbb::this_task_arena::current_thread_index();
        return idx < 0 ? 0 : static_cast<std::size_t>(idx);
    }

    /// Run `body(i)` for i in [0, n) across `workers` TBB workers.
    /// * at most one worker -> serial loop on the calling thread.
    template <class F>
    inline void for_each(std::size_t n, F &&body, std::size_t workers)
    {
        // run body over one chunk
        auto chunk = [&](const detail::Range &r)
        {
            for (std::size_t i = r.begin(); i < r.end(); ++i)
                body(i);
        };
        // single worker: whole range on the calling thread
        if (workers <= 1)
            return chunk(detail::Range(0, n));
        detail::in_arena_of(workers, [&]
                            { tbb::parallel_for(detail::Range(0, n), chunk); });
    }

    /// Run `body(i)` for i in [0, n) across `workers` TBB workers, reducing into `T` via `operator+=`.
    /// * at most one worker -> serial loop on the calling thread.
    /// @note `T` value-initializes to the identity.
    template <class T, class F>
    inline T for_reduce(std::size_t n, F &&body, std::size_t workers)
    {
        // accumulate body over one sub-range into acc
        auto chunk = [&](const detail::Range &r, T acc)
        {
            for (std::size_t i = r.begin(); i < r.end(); ++i)
                acc += body(i);
            return acc;
        };
        // single worker: whole range on the calling thread
        if (workers <= 1)
            return chunk(detail::Range(0, n), T{});
        // merge two partial sums
        auto join = [](T a, const T &b)
        { return a += b; };

        return detail::in_arena_of(
            workers, [&]
            { return tbb::parallel_reduce(detail::Range(0, n), T{}, chunk, join); });
    }

    /// Run `body(i)` for i in [0, n) across `workers`, each unit in its own arena of `inner` slots.
    ///
    /// * unit thread -> slot 0 of its arena, so current_thread_index() stays below `inner`;
    /// * at most one worker -> units in order on the calling thread;
    /// * first unit to throw -> raises `stop_flag`, a unit not yet started skips its body;
    /// * unit throwing Cancelled -> not recorded.
    /// @return Exception of the first unit to throw, with its index.
    template <class F>
    inline Failure for_units(std::size_t n,
                             F &&body,
                             std::size_t workers, std::size_t inner,
                             std::atomic<bool> &stop_flag)
    {
        Failure failure;
        // run unit i in its own arena, the first exception raising `stop`
        auto unit = [&](std::size_t i)
        {
            if (stop_flag.load(std::memory_order_relaxed))
                return;
            try
            {
                detail::in_arena(inner, [&]
                                 { body(i); });
            }
            catch (const Cancelled &)
            {
            }
            catch (...)
            {
                if (!stop_flag.exchange(true))
                    failure = {std::current_exception(), i};
            }
        };
        for_each(n, unit, workers);
        return failure;
    }

    /// Run `work()` on a new thread, calling `monitor()` on the calling thread every 250 ms.
    ///
    /// * `monitor()` throws -> raises `stop_flag`, waits for `work()`, then rethrows;
    /// * `work()` throws -> rethrown on the calling thread.
    /// @return Value of `work()`.
    template <class F, class M>
    inline auto monitored(F &&work,
                          M &&monitor,
                          std::atomic<bool> &stop_flag)
    {
        constexpr std::chrono::milliseconds period(250);
        auto result = std::async(std::launch::async, std::forward<F>(work));
        try
        {
            while (result.wait_for(period) != std::future_status::ready)
                monitor();
        }
        catch (...)
        {
            // setting shared stop_flag causes work() to interrupt
            stop_flag.store(true);
            throw;
        }
        return result.get();
    }
}
