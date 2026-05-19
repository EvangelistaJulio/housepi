/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <limits>

namespace housepi
{
    /// Split Mix 64. Used for seeding.
    /// Based on @cite vigna_pseudo_random_2015
    class SplitMix64
    {
    public:
        using result_type = std::uint64_t;
        result_type state;

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        SplitMix64() = default;
        explicit SplitMix64(result_type seed) { this->seed(seed); }

        /////////////////////////
        //////  OPERATORS  //////
        /////////////////////////

        result_type operator()() noexcept
        {
            result_type z = (state += 0x9e3779b97f4a7c15ULL);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            return z ^ (z >> 31);
        }

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        void seed(result_type seed) { state = seed; }
    };

#if __x86_64__ || __ppc64__ || __aarch64__ || _WIN64
    /// Xor Shift Rotate 256 Star Star.
    /// Based on @cite blackman_xoshiro256starstar_2018
    class XoShiRo256ss
    {
    public:
        using result_type = std::uint64_t;

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        XoShiRo256ss() = default;
        explicit XoShiRo256ss(result_type seed) { this->seed(seed); }

        /////////////////////////
        //////  OPERATORS  //////
        /////////////////////////

        result_type operator()() noexcept
        {
            const result_type r = rotl_(s_[1] * 5, 7) * 9;
            const result_type t = s_[1] << 17;

            s_[2] ^= s_[0];
            s_[3] ^= s_[1];
            s_[1] ^= s_[2];
            s_[0] ^= s_[3];

            s_[2] ^= t;
            s_[3] = rotl_(s_[3], 45);

            return r;
        }

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Seed with SplitMix64.
        void seed(result_type seed)
        {
            SplitMix64 sm(seed);
            for (result_type &x : s_)
                x = sm();
        }

        /// Jump of 2^128 steps (for parallel stream splitting).
        void jump() noexcept
        {
            static constexpr result_type J[4] = {
                0x180ec6d33cfd0abaULL,
                0xd5a61266f0c9392cULL,
                0xa9582618e03fc9aaULL,
                0x39abdc4529b1661cULL};

            result_type s0 = 0, s1 = 0, s2 = 0, s3 = 0;

            for (int i = 0; i < 4; ++i)
                for (int b = 0; b < 64; ++b)
                {
                    if (J[i] & (result_type(1) << b))
                    {
                        s0 ^= s_[0];
                        s1 ^= s_[1];
                        s2 ^= s_[2];
                        s3 ^= s_[3];
                    }
                    (*this)();
                }

            s_[0] = s0;
            s_[1] = s1;
            s_[2] = s2;
            s_[3] = s3;
        }

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        static constexpr result_type min() noexcept { return std::numeric_limits<result_type>::min(); }
        static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

    private:
        result_type s_[4]{};

        static result_type rotl_(result_type x, int k) noexcept { return (x << k) | (x >> (64 - k)); }
    };
    using rnd_gen = XoShiRo256ss;
#else

    /// Xor Shift Rotate 128 Star Star.
    /// @cite blackman_xoshiro128starstar_2018
    class XoShiRo128ss
    {
    public:
        using result_type = std::uint32_t;

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        XoShiRo128ss() = default;
        explicit XoShiRo128ss(result_type seed) { this->seed(seed); }

        /////////////////////////
        //////  OPERATORS  //////
        /////////////////////////

        result_type operator()() noexcept
        {
            const result_type r = rotl_(s_[1] * 5u, 7) * 9u;

            const result_type t = s_[1] << 9;

            s_[2] ^= s_[0];
            s_[3] ^= s_[1];
            s_[1] ^= s_[2];
            s_[0] ^= s_[3];

            s_[2] ^= t;
            s_[3] = rotl_(s_[3], 11);

            return r;
        }

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Seed with SplitMix64.
        void seed(result_type seed)
        {
            SplitMix64 sm(seed);
            std::uint64_t z = sm();
            s_[0] = static_cast<std::uint32_t>(z);
            s_[1] = static_cast<std::uint32_t>(z >> 32);
            z = sm();
            s_[2] = static_cast<std::uint32_t>(z);
            s_[3] = static_cast<std::uint32_t>(z >> 32);
        }

        /// Jump of 2^64 steps (for parallel stream splitting).
        void jump() noexcept
        {
            static constexpr result_type J[4] = {
                0x8764000bu,
                0xf542d2d3u,
                0x6fa035c3u,
                0x77f2db5bu};

            result_type s0 = 0, s1 = 0, s2 = 0, s3 = 0;

            for (int i = 0; i < 4; ++i)
                for (int b = 0; b < 32; ++b)
                {
                    if (J[i] & (result_type(1) << b))
                    {
                        s0 ^= s_[0];
                        s1 ^= s_[1];
                        s2 ^= s_[2];
                        s3 ^= s_[3];
                    }
                    (*this)();
                }

            s_[0] = s0;
            s_[1] = s1;
            s_[2] = s2;
            s_[3] = s3;
        }

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        static constexpr result_type min() noexcept { return std::numeric_limits<result_type>::min(); }
        static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

    private:
        result_type s_[4]{};

        static result_type rotl_(result_type x, int k) noexcept { return (x << k) | (x >> (32 - k)); }
    };

    using rnd_gen = XoShiRo128ss;
#endif
}
