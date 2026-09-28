/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <random>

#include <boost/random/xoshiro.hpp>

namespace housepi
{
    /// Xor Shift Rotate 256 Star Star, seeded with Split Mix 64.
    /// Based on @cite blackman_xoshiro256starstar_2018 and @cite vigna_pseudo_random_2015
    struct rnd_gen
    {
        using result_type = std::uint64_t;

        //////  Constructor  //////

        rnd_gen() = default;
        explicit rnd_gen(result_type seed) { this->seed(seed); }

        //////  Operators  //////

        /// Draw the next 64-bit integer.
        result_type operator()() noexcept { return engine_(); }

        /// Draw a standard N(0,1) deviate.
        /// Deviates come in pairs: the sampler holds a spare between calls.
        double draw_normal() { return normal_dist_(*this); }
        /// Reset the normal sampler, discarding the spare.
        void reset() { normal_dist_.reset(); }

        //////  Setter  //////

        /// Seed with SplitMix64.
        void seed(result_type seed)
        {
            reset();
            engine_.seed(seed);
        }

        /// Jump 2^128 steps ahead, splitting parallel streams.
        void jump() noexcept
        {
            reset();
            engine_.jump();
        }

        //////  Getter  //////

        /// Range of the integers drawn.
        static constexpr result_type min() noexcept { return boost::random::xoshiro256mm::min(); }
        static constexpr result_type max() noexcept { return boost::random::xoshiro256mm::max(); }

    private:
        boost::random::xoshiro256mm engine_;           /// Generator engine.
        std::normal_distribution<double> normal_dist_; /// Holds the spare of each normal pair.
    };

    //////  Standard samplers  //////

    /// @todo Improve the standard samplers.

    /// Draw a standard N(0,1) variate.
    /// The sampler lives on the generator, so seed() and jump() reset it.
    inline double standard_normal(rnd_gen &prng) { return prng.draw_normal(); }

    /// Draw a uniform [0, 1) variate, stateless.
    inline double standard_uniform(rnd_gen &prng)
    {
        std::uniform_real_distribution<double> gen(0., 1.);
        return gen(prng);
    }

    /// Draw a unit-rate exponential variate, stateless.
    inline double standard_exponential(rnd_gen &prng)
    {
        std::exponential_distribution<double> gen(1.);
        return gen(prng);
    }
}
