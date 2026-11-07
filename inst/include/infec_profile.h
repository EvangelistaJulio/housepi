/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "tools/dictionary.h"
#include "augmented.h"
#include "data/individual.h"
#include "distribution.h"

//////////////////////////////  InfecProfile  //////////////////////////////

namespace housepi
{
    namespace infprof
    {
        /// Latent instantaneous infectivity profile.
        constexpr const char *LATENT_FUNC = ".infec_prof_latent_func";
        /// Default size of the discrete-time lookup table.
        constexpr unsigned LOOKUP_SIZE = 30;

        /// Discrete-time lookup table of the generic profile.
        struct LookUpTable
        {
            //////  Setter  //////

            /// Initialize the table of `n` days from profile `ptr`.
            void init(const Distribution *ptr, size_t n);
            /// Clear the table.
            void clear() noexcept;
            /// Recover the previous table, if the generic profile is not constant.
            void recover() noexcept
            {
                if (generic_ptr_->is_constant())
                    return;
                std::swap(value_, old_);
            }
            /// Update the table, if the generic profile is not constant.
            void update()
            {
                if (generic_ptr_->is_constant())
                    return;
                std::swap(value_, old_);
                materialize_();
            }

            //////  Getter  //////

            /// Instantaneous profile on day `t`, 0 past the table.
            double inst(size_t t) const { return t < size_ ? value_.inst[t] : 0.; }
            /// Cumulative profile on day `t`, the last value past the table.
            double cumul(size_t t) const { return t < size_ ? value_.cumul[t] : value_.cumul[size_ - 1]; }

        private:
            /// Number of days in the table.
            size_t size_ = 0;
            /// Generic profile, not owned.
            const Distribution *generic_ptr_;

            /// Profile values per day.
            struct Table_
            {
                /// Instantaneous.
                std::vector<double> inst = {};
                /// Cumulative.
                std::vector<double> cumul = {};

                /// Clear the table.
                void clear() noexcept;
                /// Allocate `n` days, set to `v`.
                void assign(size_t n, double v);
            } value_, old_;

            /// Recompute every value of the table.
            void materialize_();
        };
    }

    struct Translator;

    /// Infectivity profile: generic, per-member latent, or discrete-time lookup table.
    struct InfecProfile
    {
        /// Processes the latent infectivity profile.
        friend struct Augmented;

    public:
        std::unique_ptr<Distribution> generic;
        str_table latent_definitions = {};
        size_t latent_key;

        //////  Constructor  //////

        InfecProfile();
        ~InfecProfile();

        //////  Setter  //////

        /// Build the infectivity profile distribution.
        /// @todo Honour the minimum of the profile.
        void build(Translator &translator, const str_table &input);

        /// Bind latent storage for runtime per-member profile resolution.
        void bind_augm(const Augmented *augm) noexcept { paugm_ = augm; }

        /// Finalize the generic distribution, resolving its read-surface, before use.
        /// * latent profile -> finalized with its member variable.
        void finalize()
        {
            if (generic)
                generic->finalize();
        }

        /// Store the values of a lookup-table profile.
        /// * day 0 -> infectivity 0.
        void update_cache()
        {
            switch (mode_)
            {
            case Mode::LookUp:
                lookup_.update();
                break;
            // cases Mode::Generic, Mode::Latent
            default:
                return;
            }
        }

        /// Recover the previous values of a lookup-table profile.
        void recover() noexcept
        {
            switch (mode_)
            {
            case Mode::LookUp:
                lookup_.recover();
                break;
            // cases Mode::Generic, Mode::Latent
            default:
                return;
            }
        }

        /// Request a discrete profile.
        void request_discrete() noexcept { is_required_discrete_ = true; };

        //////  Assertion  //////

        /// Throw if profile `f` transmits at generation time 0.
        /// * bound reading an undrawn value -> judged on the interval its codomain reaches;
        /// * registered profile in continuous time -> trusted.
        void assert_zero_at_origin(const Function &f) const;
        /// Throw if the generic profile transmits at generation time 0.
        void assert_valid() const
        {
            if (generic)
                assert_zero_at_origin(*generic);
        }

        //////  Getter  //////

        /// Instantaneous profile of member `indiv_ptr` at time `t`.
        double inst(const Individual *indiv_ptr, double t) const
        {
            switch (mode_)
            {
            case Mode::Latent:
                return paugm_->at(latent_key, *indiv_ptr).function()->pdf(t);
            case Mode::LookUp:
                return lookup_.inst(t);
            // case Mode::Generic
            default:
                return generic->pdf(t);
            }
        }
        /// Cumulative profile of member `indiv_ptr` at time `t`.
        double cumul(const Individual *indiv_ptr, double t) const
        {
            switch (mode_)
            {
            case Mode::Latent:
                return paugm_->at(latent_key, *indiv_ptr).function()->cdf(t);
            case Mode::LookUp:
                return generic->is_below(t) ? 0. : lookup_.cumul(t);
            // case Mode::Generic
            default:
                return generic->cdf(t);
            }
        }

        //////  Flags  //////

        /// Whether the profile is read from a latent.
        bool is_latent() const noexcept { return mode_ == Mode::Latent; }
        /// Whether a discrete profile was requested.
        bool is_required_discrete() const noexcept { return is_required_discrete_; }

    private:
        /// Source of the profile.
        enum class Mode : std::uint8_t
        {
            /// Profile read from a latent.
            Latent,
            /// Profile evaluated per call.
            Generic,
            /// Discrete-time table stored by update_cache().
            LookUp
        } mode_ = Mode::Generic;

        /// Whether a discrete profile was requested.
        bool is_required_discrete_ = false;
        /// Discrete-time lookup table.
        infprof::LookUpTable lookup_;

        /// Latent storage for the per-member profile.
        const Augmented *paugm_ = nullptr;

        /// Build the generic profile, independent of status.
        void build_generic_(Translator &translator, const str_table &input);
        /// Build the per-member profile, read from a latent function.
        void build_latent_(Translator &translator, const str_table &input);
        /// Compose a profile distribution `func` and bind its arguments.
        static std::unique_ptr<Distribution> helper_build_(Translator &translator,
                                                           const std::string &func,
                                                           bool discr);
    };
}
