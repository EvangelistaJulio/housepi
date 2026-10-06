/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
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
            void recover() noexcept;
            /// Update the table, if the generic profile is not constant.
            void update();

            //////  Getter  //////

            /// Instantaneous profile on day `t`, 0 past the table.
            double inst(size_t t) const { return t < size_ ? value_.inst[t] : 0.; }
            /// Cumulative profile on day `t`, the last value past the table.
            double cumul(size_t t) const { return t < size_ ? value_.cumul[t] : value_.cumul[size_ - 1]; }

        private:
            size_t size_ = 0;                 /// Number of days in the table.
            const Distribution *generic_ptr_; /// Generic profile, not owned.

            /// Profile values per day.
            struct Table_
            {
                std::vector<double> inst = {};  /// Instantaneous.
                std::vector<double> cumul = {}; /// Cumulative.

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
        friend struct Augmented; /// Processes the latent infectivity profile.

    public:
        std::unique_ptr<Distribution> generic; /// Generic profile.
        str_table latent_definitions = {};     /// Latent definitions.
        size_t latent_key;                     /// Latent key.

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
        /// - evaluation only: the latent profile is finalized with its variable.
        void finalize()
        {
            if (generic)
                generic->finalize();
        }

        /// Store the values of a lookup-table profile.
        /// - day 0: infectivity 0.
        void update_cache()
        {
            switch (mode_)
            {
            case Mode::LookUp:
                lookup_.update();
                break;
            default: // cases Mode::Generic, Mode::Latent
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
            default: // cases Mode::Generic, Mode::Latent
                return;
            }
        }

        /// Request a discrete profile.
        void request_discrete() noexcept { is_required_discrete_ = true; };

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
            default: // case Mode::Generic
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
            default: // case Mode::Generic
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
            Latent,  /// Profile read from a latent.
            Generic, /// Profile evaluated per call.
            LookUp   /// Discrete-time table stored by update_cache().
        } mode_ = Mode::Generic;

        bool is_required_discrete_ = false; /// Whether a discrete profile was requested.
        infprof::LookUpTable lookup_;       /// Discrete-time lookup table.

        const Augmented *paugm_ = nullptr; /// Latent storage for the per-member profile.

        /// Build the generic profile, independent of status.
        void build_generic_(Translator &translator, const str_table &input);
        /// Build the per-member profile, read from a latent function.
        void build_latent_(Translator &translator, const str_table &input);
        /// Compose a profile distribution `func` and bind its arguments.
        std::unique_ptr<Distribution> helper_build_(Translator &translator,
                                                    const std::string &func,
                                                    bool discr);
    };
}
