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
    constexpr const char *INSTANT = ".profile_inst";
    /// Default size for lookup table for discrete case
    constexpr unsigned LOOKUP_SIZE = 30;
}

struct Parameter;
class Individual;
class Distribution;

struct InfecProfile
{
    friend struct Augmented;

public:
    std::unique_ptr<Distribution> generic = nullptr;
    size_t augm_key; /// augmented profile key

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
    void update();
    void recover();
    void request_discrete() noexcept { is_discrete_requested_ = true; };

    //////////////////////
    //////  GETTER  //////
    //////////////////////

    double inst(const Individual *member_ptr, double t);
    double cumul(const Individual *member_ptr, double t);

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    bool is_augmented() const noexcept { return mode_ == Mode::Augmented; }
    bool is_lookup() const noexcept { return mode_ == Mode::LookUp; }
    bool is_discrete_requested() const noexcept { return is_discrete_requested_; }

private:
    Parameter *pParameter_;                        /// pointer to external Parameter struct
    bool is_discrete_requested_ = false;           /// wheter is requested as discrete
    size_t size_lut_ = 0;                          /// size of the LookUp Table table
    std::vector<double> lookup_inst_ = {};         /// instantaneous LookUp Table
    std::vector<double> lookup_cumul_ = {};        /// cumulative LookUp Table
    std::vector<double> lookup_buffer_inst_ = {};  /// instantaneous LookUp Table (buffer)
    std::vector<double> lookup_buffer_cumul_ = {}; /// cumulative LookUp Table (buffer)
    hsp::str_table augm_profile_definitions_ = {}; /// augmented profile definitions

    enum class Mode : uint8_t
    {
        Augmented, /// augmented infectivity profile
        Generic,   /// generic infectivity profile
        LookUp     /// generic infectivity profile
    } mode_ = Mode::Generic;

    void init_lookup_();
    void clear_lookup_();

    void build_generic_profile_(const hsp::str_table &input);
    void build_augmented_profile_(const hsp::str_table &input);
    std::unique_ptr<Distribution> helper_build_profile_(const hsp::str_table &input);
};
