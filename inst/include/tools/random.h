/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <random>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/xoshiro.hpp>

namespace housepi
{
    /// Xor Shift Rotate 256 Star Star, seeded with Split Mix 64.
    /// Based on @cite blackman_xoshiro256starstar_2018 and @cite steele_fast_2014
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

        /// Restart from the 256-bit state `state`.
        /// @note The all-zero state draws zeros only.
        void restart(const std::array<result_type, 4> &state) { engine_.seed(state.begin(), state.end()); }

        /// Jump 2^128 steps ahead, splitting parallel streams.
        void jump() noexcept { engine_.jump(); }

        //////  Getter  //////

        /// Range of the integers drawn.
        static constexpr result_type min() noexcept { return boost::random::xoshiro256mm::min(); }
        static constexpr result_type max() noexcept { return boost::random::xoshiro256mm::max(); }

    private:
        boost::random::xoshiro256mm engine_;
    };

    //////  Keyed states  //////

    namespace detail
    {
        /// SplitMix64 finalizer, a bijection of 64-bit words.
        /// Based on @cite steele_fast_2014
        inline constexpr std::uint64_t mix64(std::uint64_t z) noexcept
        {
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            return z ^ (z >> 31);
        }
    }

    /// Return the generator state of key (`seed`, `chain`, `step`, `house`), a bijection of the key.
    /// * forward pass -> word `i` mixes key word `i`, the previous word and a Weyl offset;
    /// * backward pass -> word `i` mixes itself and the next word.
    /// Every state word then reads the whole key.
    /// @note The key of the all-zero state has a seed word above 2^31, out of reach of a 32-bit seed.
    inline std::array<std::uint64_t, 4> keyed_state(std::uint64_t seed, std::uint64_t chain,
                                                    std::uint64_t step, std::uint64_t house) noexcept
    {
        constexpr std::uint64_t GAMMA = 0x9e3779b97f4a7c15ULL;
        const std::array<std::uint64_t, 4> key = {seed, chain, step, house};
        std::array<std::uint64_t, 4> state;
        std::uint64_t prev = 0;
        for (int i = 0; i < 4; ++i)
            prev = state[i] = detail::mix64(key[i] + prev + GAMMA * (i + 1));
        for (int i = 2; i >= 0; --i)
            state[i] = detail::mix64(state[i] + state[i + 1]);
        return state;
    }

    //////  Standard samplers  //////

    /// @todo Improve the standard samplers.

    /// Draw a standard N(0,1) variate with the ziggurat method, stateless.
    inline double standard_normal(RndGen &prng) { return boost::random::normal_distribution<double>()(prng); }

    /// Draw a uniform [0, 1) variate, stateless.
    inline double standard_uniform(RndGen &prng) { return std::uniform_real_distribution<double>()(prng); }

    /// Draw a unit-rate exponential variate, stateless.
    inline double standard_exponential(RndGen &prng) { return std::exponential_distribution<double>()(prng); }

    /// Draw a unit-scale gamma variate of shape `shape`, stateless.
    /// Based on @cite marsaglia_simple_2000
    /// * shape at least 1 -> squeeze and rejection on a cubed normal;
    /// * shape below 1 -> variate of shape + 1, times U^(1 / shape).
    inline double standard_gamma(RndGen &prng, double shape)
    {
        if (shape < 1.)
            return standard_gamma(prng, shape + 1.) * std::pow(standard_uniform(prng), 1. / shape);

        const double d = shape - 1. / 3., c = 1. / std::sqrt(9. * d);
        double z, z2, v, u;
        do
        {
            do
            {
                z = standard_normal(prng);
                v = 1. + c * z;
            } while (v <= 0.);
            v = v * v * v;
            z2 = z * z;
            u = standard_uniform(prng);
        } while (u >= 1. - 0.0331 * z2 * z2 && std::log(u) >= 0.5 * z2 + d * (1. - v + std::log(v)));
        return d * v;
    }
}
