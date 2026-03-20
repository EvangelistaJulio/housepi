/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "tools/types.h"

#include "distribution.h"
#include "individual.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** InfecProfile : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp
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

                void clear() noexcept;
                void assign(size_t n, double v);
            } value_, old_;

            /// Recompute all values of the table.
            void materialize_();
        };
    }

    struct Parameter;

    struct InfecProfile
    {
        friend struct Augmented;

    public:
        std::unique_ptr<Distribution> generic; /// generic profile
        str_table augm_definitions = {};       /// augmented definitions
        size_t augm_key;                       /// augmented key

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        InfecProfile(Parameter *paramptr);
        ~InfecProfile();

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        void build(const str_table &input);
        void init_template(str_table input);
        void init_template(dict<std::string> input);

        /// If infectivity profile is a lookup table, store its values.
        /// NOTE: infectivity is set to 0 on day 0.
        void update()
        {
            switch (mode_)
            {
            case Mode::Generic:
            case Mode::Augmented:
                return;
            case Mode::LookUp:
                lookup_.update();
                break;
            default:
                __builtin_unreachable();
            }
        }

        /// If infectivity profile is a lookup table, recover old values.
        void recover()
        {
            switch (mode_)
            {
            case Mode::Generic:
            case Mode::Augmented:
                return;
            case Mode::LookUp:
                lookup_.recover();
                break;
            default:
                __builtin_unreachable();
            }
        }

        void request_discrete() noexcept { is_discrete_requested_ = true; };

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        double inst(const Individual *member_ptr, double t) const
        {
            switch (mode_)
            {
            case Mode::Augmented:
                return member_ptr->get_augm_variable(augm_key).function()->pdf(t);
            case Mode::Generic:
                return generic->pdf(t);
            case Mode::LookUp:
                return lookup_.inst(t);
            }
            __builtin_unreachable();
        }
        double cumul(const Individual *member_ptr, double t) const
        {
            switch (mode_)
            {
            case Mode::Augmented:
                return member_ptr->get_augm_variable(augm_key).function()->cdf(t);
            case Mode::Generic:
                return generic->cdf(t);
            case Mode::LookUp:
                return generic->is_below(t) ? 0. : lookup_.cumul(t);
            }
            __builtin_unreachable();
        }

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        bool is_augmented() const noexcept { return mode_ == Mode::Augmented; }
        bool is_lookup() const noexcept { return mode_ == Mode::LookUp; }
        bool is_discrete_requested() const noexcept { return is_discrete_requested_; }

    private:
        enum class Mode : std::uint8_t
        {
            Augmented, /// augmented infectivity profile
            Generic,   /// generic infectivity profile
            LookUp     /// generic infectivity profile
        } mode_ = Mode::Generic;

        bool is_discrete_requested_ = false; /// wheter is requested as discrete
        Parameter *pParameter_;              /// pointer to external Parameter struct
        infprof::LookUpTable lookup_;        /// discrete time lookup table

        /// Build for generic infectivity profile (no status dependence).
        void build_generic_(const str_table &input);
        /// Build for individual augmented function infectivity profile.
        void build_augmented_(const str_table &input);
        /// Helper to build profile.
        std::unique_ptr<Distribution> helper_build_(const str_table &input);
    };
}
