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

#include "tools/types.h"

#include "augmented.h"
#include "distribution.h"
#include "individual.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** InfecProfile : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    namespace infprof
    {
        /// Augmented instantaneous infectivity profile.
        constexpr const char *AUGM_FUNC = ".infec_prof_augm_func";
        /// Default size for lookup table for discrete case
        constexpr unsigned LOOKUP_SIZE = 30;

        /// Helper module to deal with the case of discrete time lookup table.
        struct LookUpTable
        {
            //////////////////////
            //////  SETTER  //////
            //////////////////////

            /// Initialize lookuup table.
            void init(const Distribution *ptr, size_t n);
            /// Clear infectivity profile lookup table.
            void clear() noexcept;
            /// Recover old state if generic is not constant.
            void recover() noexcept;
            /// Update lookup table if generic is not constant.
            void update();

            //////////////////////
            //////  GETTER  //////
            //////////////////////

            double inst(size_t t) const { return t < size_ ? value_.inst[t] : 0.; }
            double cumul(size_t t) const { return t < size_ ? value_.cumul[t] : value_.cumul[size_ - 1]; }

        private:
            size_t size_ = 0;                 /// size of look up table
            const Distribution *generic_ptr_; /// non-owning pointer to generic

            struct Table_
            {
                std::vector<double> inst = {};  /// instantaneous
                std::vector<double> cumul = {}; /// cumulative

                /// Clear lookup table.
                void clear() noexcept;
                /// Allocate lookup table.
                void assign(size_t n, double v);
            } value_, old_;

            /// Recompute all values of the table.
            void materialize_();
        };
    }

    struct Translator;

    struct InfecProfile
    {
        friend struct Augmented; /// process augmented infectivity profile

    public:
        std::unique_ptr<Distribution> generic; /// generic profile
        str_table augm_definitions = {};       /// augmented definitions
        size_t augm_key;                       /// augmented key

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        InfecProfile();
        ~InfecProfile();

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Build infectivity profile distribution.
        /// NOTE: at the moment I am ignoring minimum
        void build(Translator &translator, const str_table &input);

        /// Bind augmented storage for runtime per-member profile resolution.
        void bind_augm(const Augmented *augm) noexcept { paugm_ = augm; }

        /// Finalize the generic distribution (resolve its read-surface) before use.
        /// NOTE: evaluation-only; the augmented profile is finalized with its variable.
        void finalize()
        {
            if (generic)
                generic->finalize();
        }

        /// If infectivity profile is a lookup table, store its values.
        /// NOTE: infectivity is set to 0 on day 0.
        void update_cache()
        {
            switch (mode_)
            {
            case Mode::LookUp:
                lookup_.update();
                break;
            default: // cases Mode::Generic, Mode::Augmented
                return;
            }
        }

        /// If infectivity profile is a lookup table, recover old values.
        void recover() noexcept
        {
            switch (mode_)
            {
            case Mode::LookUp:
                lookup_.recover();
                break;
            default: // cases Mode::Generic, Mode::Augmented
                return;
            }
        }

        void request_discrete() noexcept { is_required_discrete_ = true; };

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        double inst(const Individual *indiv_ptr, double t) const
        {
            switch (mode_)
            {
            case Mode::Augmented:
                return paugm_->at(augm_key, *indiv_ptr).function()->pdf(t);
            case Mode::LookUp:
                return lookup_.inst(t);
            default: // case Mode::Generic
                return generic->pdf(t);
            }
        }
        double cumul(const Individual *indiv_ptr, double t) const
        {
            switch (mode_)
            {
            case Mode::Augmented:
                return paugm_->at(augm_key, *indiv_ptr).function()->cdf(t);
            case Mode::LookUp:
                return generic->is_below(t) ? 0. : lookup_.cumul(t);
            default: // case Mode::Generic
                return generic->cdf(t);
            }
        }

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        bool is_augmented() const noexcept { return mode_ == Mode::Augmented; }
        bool is_lookup() const noexcept { return mode_ == Mode::LookUp; }
        bool is_required_discrete() const noexcept { return is_required_discrete_; }

    private:
        enum class Mode : std::uint8_t
        {
            Augmented, /// augmented infectivity profile
            Generic,   /// generic infectivity profile
            LookUp     /// generic infectivity profile
        } mode_ = Mode::Generic;

        bool is_required_discrete_ = false; /// wheter is requested as discrete
        infprof::LookUpTable lookup_;       /// discrete time lookup table

        const Augmented *paugm_ = nullptr; /// augmented storage for per-member profile

        /// Build for generic infectivity profile (no status dependence).
        void build_generic_(Translator &translator, const str_table &input);
        /// Build for individual augmented function infectivity profile.
        void build_augmented_(Translator &translator, const str_table &input);
        /// Helper to compose a profile distribution and bind its arguments.
        std::unique_ptr<Distribution> helper_build_(Translator &translator,
                                                    const std::string &func,
                                                    bool discr);
    };
}
