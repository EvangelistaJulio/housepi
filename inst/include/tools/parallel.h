/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <exception>
#include <future>
#include <utility>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

namespace housepi::parallel
{
    /// First error of a run of units.
    struct Failure
    {
        /// Exception that stopped the run, null where every unit completed.
        std::exception_ptr error;
        /// Unit that raised it.
        std::size_t index = 0;
    };

    /// Thrown by a unit that stops because another unit failed.
    struct Cancelled
    {
    };

    namespace detail
    {
        using Range = tbb::blocked_range<std::size_t>;

        /// Set x87 precision of the calling thread to 64 bits, as on R's threads.
        /// * Windows x86-64 -> TBB workers start at 53 bits, MinGW's libm computes on the x87;
        /// * elsewhere -> no-op.
        inline void x87_extended_precision() noexcept
        {
#if defined(_WIN32) && defined(__x86_64__)
            unsigned short cw;
            __asm__ volatile("fnstcw %0" : "=m"(cw));
            if ((cw & 0x300) != 0x300)
            {
                cw |= 0x300;
                __asm__ volatile("fldcw %0" : : "m"(cw));
            }
#endif
        }

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

        /// Sum `body(i)` for i in [begin, end), left to right.
        template <class T, class F>
        inline T sum_range(F &body, std::size_t begin, std::size_t end)
        {
            T acc{};
            for (std::size_t i = begin; i < end; ++i)
                acc += body(i);
            return acc;
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
    /// * at most one worker -> serial loop on the calling thread;
    /// * any chunk -> 64-bit x87 precision on the thread running it.
    template <class F>
    inline void for_each(std::size_t n, F &&body, std::size_t workers)
    {
        // run body over one chunk
        auto chunk = [&](const detail::Range &r)
        {
            detail::x87_extended_precision();
            for (std::size_t i = r.begin(); i < r.end(); ++i)
                body(i);
        };
        // single worker: whole range on the calling thread
        if (workers <= 1)
            return chunk(detail::Range(0, n));
        detail::in_arena_of(workers, [&]
                            { tbb::parallel_for(detail::Range(0, n), chunk); });
    }

    /// Maximum number of blocks a reduction cuts its range into.
    constexpr std::size_t REDUCE_BLOCKS = 128;

    /// Run `body(i)` for i in [0, n) across `workers` TBB workers, reducing into `T` via `operator+=`.
    ///
    /// Same result at any `workers`, 1 included:
    /// * range -> `min(n, REDUCE_BLOCKS)` contiguous blocks, each summed left to right;
    /// * block sums -> added in block order on the calling thread;
    /// * at most one worker -> blocks in order on the calling thread;
    /// * any chunk -> 64-bit x87 precision on the thread running it.
    /// @note `T` value-initializes to the identity.
    template <class T, class F>
    inline T for_reduce(std::size_t n, F &&body, std::size_t workers)
    {
        const std::size_t n_blocks = std::min(n, REDUCE_BLOCKS);
        // first index of block b, the first n % n_blocks blocks one term longer
        const std::size_t q = n / n_blocks, r = n % n_blocks;
        auto start = [&](std::size_t b)
        { return b * q + std::min(b, r); };
        T total{};
        // single worker: blocks in order on the calling thread
        if (workers <= 1)
        {
            detail::x87_extended_precision();
            for (std::size_t b = 0; b < n_blocks; ++b)
                total += detail::sum_range<T>(body, start(b), start(b + 1));
            return total;
        }
        std::array<T, REDUCE_BLOCKS> sums;
        for_each(n_blocks, [&](std::size_t b)
                 { sums[b] = detail::sum_range<T>(body, start(b), start(b + 1)); }, workers);
        for (std::size_t b = 0; b < n_blocks; ++b)
            total += sums[b];
        return total;
    }

    /// Run `body(i)` for i in [0, n) across `workers`, each unit in its own arena of `inner` slots.
    ///
    /// * units -> dispatched in an arena of `workers` slots of their own, never the caller's;
    /// * unit thread -> slot 0 of its arena, so current_thread_index() stays below `inner`;
    /// * unit -> 64-bit x87 precision on the thread running it;
    /// * at most one worker -> units in order on the calling thread;
    /// * unit throwing -> raises `stop_flag`, a unit not yet started skips its body;
    /// * first unit to throw in this call -> recorded, whoever else raised `stop_flag`;
    /// * unit throwing Cancelled -> not recorded.
    /// @return Exception of the first unit to throw, with its index.
    template <class F>
    inline Failure for_units(std::size_t n,
                             F &&body,
                             std::size_t workers, std::size_t inner,
                             std::atomic<bool> &stop_flag)
    {
        Failure failure;
        std::atomic<bool> is_failed{false};
        // run unit i in its own arena, any exception raising `stop_flag`
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
                stop_flag.store(true);
                if (!is_failed.exchange(true))
                    failure = {std::current_exception(), i};
            }
        };
        // an arena of their own, whatever the size of the caller's
        detail::in_arena(workers, [&]
                         { for_each(n, unit, workers); });
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
        constexpr std::chrono::milliseconds PERIOD(250);
        auto result = std::async(std::launch::async, std::forward<F>(work));
        try
        {
            while (result.wait_for(PERIOD) != std::future_status::ready)
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
