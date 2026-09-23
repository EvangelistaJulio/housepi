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

        /////////////////////
        ///  CONSTRUCTOR  ///
        /////////////////////

        rnd_gen() = default;
        explicit rnd_gen(result_type seed) { this->seed(seed); }

        ///////////////////
        ///  OPERATORS  ///
        ///////////////////

        result_type operator()() noexcept { return engine_(); }

        /// Draw a standard N(0,1) deviate.
        /// NOTE: deviates come in pairs, sampler holds a spare between calls.
        double draw_normal() { return normal_dist_(*this); }
        ///  Reset the normal distribution sampler (discard spare).
        void reset() { normal_dist_.reset(); }

        ////////////////
        ///  SETTER  ///
        ////////////////

        /// Seed with SplitMix64.
        void seed(result_type seed)
        {
            reset();
            engine_.seed(seed);
        }

        /// Jump of 2^128 steps (for parallel stream splitting).
        void jump() noexcept
        {
            reset();
            engine_.jump();
        }

        ////////////////
        ///  GETTER  ///
        ////////////////

        static constexpr result_type min() noexcept { return boost::random::xoshiro256mm::min(); }
        static constexpr result_type max() noexcept { return boost::random::xoshiro256mm::max(); }

    private:
        boost::random::xoshiro256mm engine_;
        std::normal_distribution<double> normal_dist_; /// Holds the spare of each normal pair.
    };

    ///////////////////////////
    ///  STANDARD SAMPLERS  ///
    ///////////////////////////

    /// UPDATE: these samplers can be improved.

    /// Standard N(0,1) random variable.
    /// NOTE: sampler lives on the generator in order to allow for reset.
    inline double standard_normal(rnd_gen &prng) { return prng.draw_normal(); }

    /// Uniform [0,1) random variable. Stateless.
    inline double standard_uniform(rnd_gen &prng)
    {
        std::uniform_real_distribution<double> gen(0., 1.);
        return gen(prng);
    }

    /// Exponential random variable with rate 1. Stateless.
    inline double standard_exponential(rnd_gen &prng)
    {
        std::exponential_distribution<double> gen(1.);
        return gen(prng);
    }
}
