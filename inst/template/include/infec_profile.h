/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__INFEC_PROFILE__H
#define HOUSEPI__INFEC_PROFILE__H

#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  *** InfecProfile : struct ***  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

struct Parameter;
class Individual;
class Distribution;

struct InfecProfile
{
public:
    std::unique_ptr<Distribution> generic_profile;

    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////  InfecProfile : constructor  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    InfecProfile(std::shared_ptr<Parameter> paramptr);
    ~InfecProfile();

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  InfecProfile : setter  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    void build(const u_dict<str_vec> &deftable);
    void init_template(u_dict<str_vec> deftable);
    void update();
    void recover();
    void request_discrete() { is_discrete_requested_ = true; };
    void transfer_augm_definitions(u_dict<str_vec> &target) const;

    /////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////  InfecProfile : getter  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    bool is_lookup() const { return inst_ == &InfecProfile::do_inst_lookup_; }
    bool is_augmented() const { return inst_ == &InfecProfile::do_inst_augm_; }
    bool is_discrete_requested() const { return is_discrete_requested_; }

    double inst(const Individual *member_ptr, double t) { return (this->*inst_)(member_ptr, t); }
    double cumul(const Individual *member_ptr, double t) { return (this->*cumul_)(member_ptr, t); }

private:
    std::shared_ptr<Parameter> pParameter_;         /// pointer to external Parameter struct
    bool is_discrete_requested_ = false;            /// wheter is requested as discrete
    size_t size_lut_ = 0;                           /// size of the LookUp Table table
    dbl_vec lookup_inst_ = {};                      /// instantaneous LookUp Table
    dbl_vec lookup_cumul_ = {};                     /// cumulative LookUp Table
    dbl_vec lookup_buffer_inst_ = {};               /// instantaneous LookUp Table (buffer)
    dbl_vec lookup_buffer_cumul_ = {};              /// cumulative LookUp Table (buffer)
    u_dict<str_vec> augm_profile_definitions_ = {}; /// augmented profile definitions

    void init_lookup_();
    void clear_lookup_();

    void build_generic_profile_(const u_dict<str_vec> &input);
    void build_augmented_profile_(const u_dict<str_vec> &input);
    std::unique_ptr<Distribution> helper_build_profile_(const u_dict<std::string> &input);

    //===== Functional methods for getter methods =====//
    double (InfecProfile::*inst_)(const Individual *, double) = &InfecProfile::do_inst_default_;
    double do_inst_default_(const Individual *, double t);
    double do_inst_lookup_(const Individual *, double t);
    double do_inst_augm_(const Individual *i, double t);

    double (InfecProfile::*cumul_)(const Individual *, double) = &InfecProfile::do_cumul_default_;
    double do_cumul_default_(const Individual *, double t);
    double do_cumul_lookup_(const Individual *, double t);
    double do_cumul_augm_(const Individual *i, double t);
};

#endif // HOUSEPI__INFEC_PROFILE__H
