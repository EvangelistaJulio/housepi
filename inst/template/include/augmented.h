/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <vector>

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
    ~Augmented();

    std::unordered_set<size_t> list_draws = {};
    std::unordered_set<size_t> list_export = {};
    std::unordered_set<size_t> user_discrete_set = {};

    /// Status to definition row.
    std::unordered_map<int, std::vector<size_t>> status_to_def_row = {};
    /// Status to key of augmented variable to be sampled.
    std::unordered_map<int, std::vector<size_t>> status_to_sample = {};
    /// Status to internal index.
    std::unordered_map<int, size_t> status_to_index = {};

    // optional augmented posterior sample
    std::vector<const Individual *> smpl_individuals = {};
    std::vector<std::vector<std::vector<double>>> smpl_values = {};

    Augmented(std::shared_ptr<DataWide> dataptr,
              std::shared_ptr<Parameter> paramptr,
              std::shared_ptr<InfecProfile> infprofptr);

    void build(const hsp::str_table &input);
    void init();
    void init_individual_sample(std::vector<std::string> &id_draws);

    void user_reset(hsp::rnd_gen &prng);
    void user_reset(hsp::rnd_gen &prng, Individual &member);

    void clear_all(Individual &member);
    void init_missing(Individual &member, bool force = false);

    void init_draws(unsigned long sample_size);
    void store_draws(unsigned long smpl_count);

    /// Flag if there are non-constant augmented definitions.
    bool is_any_sample() const { return !status_to_sample.empty(); }
    /// Flag if there are augmented definitions with a model.
    bool is_any_model() const { return is_any_model_; }
    /// Flag if augmented draws are stored.
    bool is_any_draw() const { return !smpl_individuals.empty(); }

    double materialize_augm_model_log_lik(Individual &member);
    double materialize_augm_model_logratio(Individual &member, size_t a);

private:
    std::shared_ptr<DataWide> pDataWide_;         /// pointer to external DataWide struct
    std::shared_ptr<Parameter> pParameter_;       /// pointer to external Parameter struct
    std::shared_ptr<InfecProfile> pInfecProfile_; /// pointer to external InfecProfile struct
    hsp::str_table deftable_ = {};                /// store definition table


    struct ModelBlock_
    {
        std::vector<size_t> models;
        std::unordered_map<size_t, std::vector<size_t>> depen;
    };
    std::vector<ModelBlock_> model_blocks_;
    bool is_any_model_ = false; /// store flag if one augmented model is provided

    void helper_build_definition_by_status_();
    void helper_build_definition_order_();
    void set_definition_dependency_order_for_(int status);

    std::vector<size_t> get_unsorted_definition_rows_(int status, bool is_pre_process) const;
    hsp::str_table get_ordered_models_(int status) const;

    void set_sample_map_for_(int status, size_t k);
    void set_function_or_identity_for_(size_t k);
    void init_model_list();

    void helper_init_attributes_(Individual &member, size_t k);
    std::vector<size_t> init_get_definition_order_(Individual &member);
    void init_model_attributes_(Individual &member, size_t k);
    void init_func_attributes_(Individual &member, size_t k);
};
