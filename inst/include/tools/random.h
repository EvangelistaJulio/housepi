/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>
#include <random>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/xoshiro.hpp>

namespace housepi
{
    /// Xor Shift Rotate 256 Star Star, seeded with Split Mix 64.
    /// Based on @cite blackman_xoshiro256starstar_2018 and @cite vigna_pseudo_random_2015
    struct RndGen
    {
        using result_type = std::uint64_t;

        //////  Constructor  //////

        RndGen() = default;
        explicit RndGen(result_type seed) { this->seed(seed); }

        //////  Operators  //////

        /// Draw the next 64-bit integer.
        result_type operator()() noexcept { return engine_(); }

        //////  Setter  //////

        /// Seed with SplitMix64.
        void seed(result_type seed) { engine_.seed(seed); }

        /// Jump 2^128 steps ahead, splitting parallel streams.
        void jump() noexcept { engine_.jump(); }

        //////  Getter  //////

        /// Range of the integers drawn.
        static constexpr result_type min() noexcept { return boost::random::xoshiro256mm::min(); }
        static constexpr result_type max() noexcept { return boost::random::xoshiro256mm::max(); }

    private:
        boost::random::xoshiro256mm engine_; /// Generator engine.
    };

    //////  Standard samplers  //////

    /// @todo Improve the standard samplers.

    /// Draw a standard N(0,1) variate with the ziggurat method, stateless.
    inline double standard_normal(RndGen &prng) { return boost::random::normal_distribution<double>()(prng); }

    /// Draw a uniform [0, 1) variate, stateless.
    inline double standard_uniform(RndGen &prng) { return std::uniform_real_distribution<double>()(prng); }

    /// Draw a unit-rate exponential variate, stateless.
    inline double standard_exponential(RndGen &prng) { return std::exponential_distribution<double>()(prng); }
}
