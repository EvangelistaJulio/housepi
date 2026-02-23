/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** InfecProfile : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace hsp::infprof
{
    /// Augmented instantaneous infectivity profile.
    constexpr const char *AUGM_FUNC = ".infec_prof_augm_func";
    /// Default size for lookup table for discrete case
    constexpr unsigned LOOKUP_SIZE = 30;
}

struct Parameter;
class Distribution;
class AugmContainer;

struct InfecProfile
{
    friend struct Augmented;

public:
    std::unique_ptr<Distribution> generic; /// generic profile
    size_t augm_key;                       /// augmented key
    hsp::str_table augm_definitions = {};  /// augmented definitions

    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    InfecProfile(Parameter *paramptr);
    ~InfecProfile();

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    void build(const hsp::str_table &input);
    void init_template(hsp::str_table input);
    void init_template(hsp::dict<std::string> input);

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

    /// If infectivity profile is a lookup table, recover its buffered values.
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

    double inst(const AugmContainer *member_ptr, double t) const;
    double cumul(const AugmContainer *member_ptr, double t) const;

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    bool is_augmented() const noexcept { return mode_ == Mode::Augmented; }
    bool is_lookup() const noexcept { return mode_ == Mode::LookUp; }
    bool is_discrete_requested() const noexcept { return is_discrete_requested_; }

private:
    enum class Mode : uint8_t
    {
        Augmented, /// augmented infectivity profile
        Generic,   /// generic infectivity profile
        LookUp     /// generic infectivity profile
    } mode_ = Mode::Generic;

    Parameter *pParameter_;              /// pointer to external Parameter struct
    bool is_discrete_requested_ = false; /// wheter is requested as discrete

    /// Build for generic infectivity profile.
    /// NOTE: no status dependence allowed here.
    void build_generic_(const hsp::str_table &input);
    /// Build for individual augmented function infectivity profile.
    void build_augmented_(const hsp::str_table &input);
    /// Helper to build profile.
    std::unique_ptr<Distribution> helper_build_(const hsp::str_table &input);

    /// Helper module to deal with the case of discrete time lookup table.
    struct LookUpTable
    {
        /// Initialize lookuup table.
        void init(const Distribution *ptr, size_t n)
        {
            generic_ptr_ = ptr;
            size_ = n;
            inst_.assign(n, 0);
            cumul_.assign(n, 0);
            buffer_inst_.assign(n, 0);
            buffer_cumul_.assign(n, 0);
            materialize_();
        }
        /// Clear infectivity profile lookup table.
        void clear() noexcept
        {
            inst_.clear();
            cumul_.clear();
            buffer_inst_.clear();
            buffer_cumul_.clear();
        }
        /// Recover buffered state if generic is not constant.
        void recover() noexcept
        {
            if (!is_constant())
                swap_buffer_();
        }
        /// Update lookup table if generic is not constant.
        void update()
        {
            if (!is_constant())
            {
                swap_buffer_();
                materialize_();
            }
        }

        double inst(size_t t) const { return t < size_ ? inst_[t] : 0.; }
        double cumul(size_t t) const { return t < size_ ? cumul_[t] : cumul_[size_ - 1]; }

        bool is_constant() const noexcept;

    private:
        size_t size_ = 0;                       /// size of look up table
        std::vector<double> inst_ = {};         /// instantaneous
        std::vector<double> cumul_ = {};        /// cumulative
        std::vector<double> buffer_inst_ = {};  /// instantaneous buffer
        std::vector<double> buffer_cumul_ = {}; /// cumulativ buffer
        const Distribution *generic_ptr_;       /// non-owning pointer to generic

        /// Swap current lookup table with buffer.
        void swap_buffer_() noexcept
        {
            inst_.swap(buffer_inst_);
            cumul_.swap(buffer_cumul_);
        }
        /// Recompute all values of the table.
        void materialize_();
    } lookup_;
};
