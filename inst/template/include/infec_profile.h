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

    double inst(const Individual *member_ptr, double t) { return (this->*inst_)(member_ptr, t); }
    double cumul(const Individual *member_ptr, double t) { return (this->*cumul_)(member_ptr, t); }

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    bool is_lookup() const noexcept { return inst_ == &InfecProfile::do_inst_lookup_; }
    bool is_augmented() const noexcept { return inst_ == &InfecProfile::do_inst_augm_; }
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

    void init_lookup_();
    void clear_lookup_();

    void build_generic_profile_(const hsp::str_table &input);
    void build_augmented_profile_(const hsp::str_table &input);
    std::unique_ptr<Distribution> helper_build_profile_(const hsp::str_table &input);

    /// Functional method
    double (InfecProfile::*inst_)(const Individual *, double) = &InfecProfile::do_inst_default_;
    /// Default method for instantaneous infectivity profile.
    double do_inst_default_(const Individual *, double t);
    /// Lookup method for discrete instantaneous infectivity profile.
    double do_inst_lookup_(const Individual *, double t);
    /// Augmented method for instantaneous infectivity profile.
    double do_inst_augm_(const Individual *i, double t);

    double (InfecProfile::*cumul_)(const Individual *, double) = &InfecProfile::do_cumul_default_;
    /// Default method for cumulative infectivity profile.
    double do_cumul_default_(const Individual *, double t);
    /// Lookup method for discrete cumulative infectivity profile.
    double do_cumul_lookup_(const Individual *, double t);
    /// Augmented method for cumulative infectivity profile.
    double do_cumul_augm_(const Individual *i, double t);
};
