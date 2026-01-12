/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__AUGMENTED__H
#define HOUSEPI__AUGMENTED__H

#include <memory>

#include "tools/random.h"
#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Augmented : struct ***  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

struct DataWide;
struct Parameter;
struct InfecProfile;
class Individual;

struct Augmented
{
    str_vec list_draws = {};
    str_vec list_full = {};
    str_vec list_export = {};
    str_vec user_list_discrete = {};

    u_dict<std::vector<size_t>> status_to_def_row = {}; /// definition rows by status
    std::unordered_map<int, u_dict<str_vec>> status_map_model = {};
    std::unordered_map<int, str_vec> status_map_sample = {};

    // optional augmented posterior sample
    std::vector<const Individual *> smpl_individuals = {};
    std::vector<u_dict<dbl_vec>> smpl_values = {};

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  Augmented : constructor  /////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    Augmented(std::shared_ptr<DataWide> dataptr,
              std::shared_ptr<Parameter> paramptr,
              std::shared_ptr<InfecProfile> infprofptr);
    ~Augmented();

    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : setter  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    void build(const u_dict<str_vec> &deftable);
    void init();
    void init_individual_sample(str_vec &id_draws);

    void user_reset(rnd_gen &prng);
    void user_reset(rnd_gen &prng, Individual &member);

    void clear_all(Individual &member);
    void init_all(Individual &member);
    void init_missing(Individual &member);
    void init_model_logratios(Individual &member);

    void init_draws(unsigned long sample_size);
    void store_draws(unsigned long smpl_count);

    /////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  Augmented : getter  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    /// Flag if there are any augmented definitions.
    bool is_any() const { return !list_full.empty(); }
    /// Flag if there are non-constant augmented definitions.
    bool is_any_sample() const { return !status_map_sample.empty(); }
    /// Flag if there are augmented definitions with a model.
    bool is_any_model() const { return !status_map_model.empty(); }
    /// Flag if augmented draws are stored.
    bool is_any_draw() const { return !smpl_individuals.empty(); }

private:
    std::shared_ptr<DataWide> pDataWide_;         /// pointer to external DataWide struct
    std::shared_ptr<Parameter> pParameter_;       /// pointer to external Parameter struct
    std::shared_ptr<InfecProfile> pInfecProfile_; /// pointer to external InfecProfile struct
    u_dict<str_vec> deftable_ = {};               /// store definition table

    void helper_assert_no_reserved_();
    void helper_build_definition_by_status_();
    void helper_build_definition_order_();
    void set_definition_dependency_order_for_(const std::string &s);
    std::vector<size_t> get_unsorted_definition_rows_(const std::string &s, bool is_pre_process) const;
    void set_sample_map_for_(const std::string &s, size_t k);
    void set_function_or_identity_for_(size_t k);
    void helper_build_model_status_maps_();

    void helper_init_attributes_(Individual &member, size_t k);
    std::vector<size_t> init_get_definition_order_(Individual &member);
    void init_model_attributes_(Individual &member, size_t k);
    void init_func_attributes_(Individual &member, size_t k);
};

#endif // HOUSEPI__AUGMENTED__H
