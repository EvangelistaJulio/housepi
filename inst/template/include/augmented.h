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

namespace hsp::key
{
    class Registry;
}

struct DataWide;
struct Parameter;
struct InfecProfile;
class AugmContainer;
class Individual;

struct Augmented
{
    /// UPATE: once read_attribute goes back to MCMCImpl, this can become private
    std::unique_ptr<hsp::key::Registry> pRegistry; /// pointer to registry

    std::unordered_set<size_t> list_draws = {};
    std::unordered_set<size_t> list_export = {};

    /// Status to key of augmented variable to be sampled.
    std::unordered_map<int, std::vector<size_t>> status_to_sample = {};

    // optional augmented posterior sample
    std::vector<const AugmContainer *> smpl_individuals = {};
    std::vector<std::vector<std::vector<double>>> smpl_values = {};

    ///////////////////////////
    //////  CONSTRUCTOR  //////
    ///////////////////////////

    Augmented(const std::vector<std::string> &keys,
              DataWide *dataptr,
              Parameter *paramptr,
              InfecProfile *infprofptr);
    ~Augmented();

    //////////////////////
    //////  SETTER  //////
    //////////////////////

    /// Assign augmented variables to all individuals.
    void build(const hsp::str_table &input);

    /// Assign augmented variables only if missing.
    /// If force is true, re-initialize all augmented variables.
    void init_missing(Individual &member, bool force = false);
    /// Initialize augmented dependecies for evaluations.
    /// NOTE: not to be called before simulation which builds on demand.
    void init();

    /// Request discrete augmented definition.
    void request_discrete(size_t a) { user_discrete_set_.emplace(a); }

    /// Reset augmented initial values for given individual.
    void user_reset(hsp::rnd_gen &prng, AugmContainer &member);
    /// Initialize augmented initial values for all individuals.
    void user_reset(hsp::rnd_gen &prng);

    /// Set individual id for augmented draws follow-up.
    void init_individual_sample(std::vector<std::string> &id_draws);
    /// Initialize sample for individual augmented variable followup.
    void init_draws(unsigned long sample_size);
    /// Store sample for individual augmented variable followup.
    void store_draws(unsigned long smpl_count);

    //////////////////////
    //////  GETTER  //////
    //////////////////////

    /// Compute and store log-likelihood associated to augmented variable models.
    /// NOTE: to be called only if any model is defined.
    double materialize_augm_model_log_lik(AugmContainer &member);
    /// Change in log-likelihood associated to given augmented.
    double materialize_augm_model_logratio(AugmContainer &member, size_t a);

    /////////////////////
    //////  FLAGS  //////
    /////////////////////

    /// Flag if there are non-constant augmented definitions.
    bool is_any_sample() const noexcept { return is_any_sample_; }
    /// Flag if there are augmented definitions with a model.
    bool is_any_model() const noexcept { return is_any_model_; }
    /// Flag if augmented draws are stored.
    bool is_any_draw() const noexcept { return !smpl_individuals.empty(); }

private:
    DataWide *pDataWide_;          /// pointer to external DataWide struct
    Parameter *pParameter_;        /// pointer to external Parameter struct
    InfecProfile *pInfecProfile_;  /// pointer to external InfecProfile struct
    hsp::str_table deftable_ = {}; /// store definition table

    bool is_any_sample_ = false; /// whether at least one augmented model is provided
    bool is_any_model_ = false;  /// whether at least one augmented model is provided

    /// User defined discrete parameters
    std::unordered_set<size_t> user_discrete_set_ = {};

    /// List of all infected statuses in the dataset
    std::unordered_set<int> all_infected_statuses_ = {};
    /// Status to definition row.
    std::unordered_map<int, std::vector<size_t>> status_to_def_row_ = {};
    /// Status to internal index.
    std::unordered_map<int, size_t> status_to_index_ = {};

    /// Auxiliary block for status-specific models.
    struct ModelBlock_
    {
        std::vector<size_t> models = {};             /// list of augmented with models
        std::vector<std::vector<size_t>> depen = {}; /// list of respective dependecies
    };
    std::vector<ModelBlock_> model_blocks_;

    /// Eventually append infectivity profile definitions.
    void helper_build_append_infec_profile_();

    /// Assert that augmented status/generic definitions are non degenerate.
    /// Also, assert that status-specific definitions are for infected status only.
    /// By convention, generic defintion are associated to missing status.
    void helper_build_definition_by_status_();

    /// Emplace all augmented definitions.
    void helper_build_emplace_all_();

    /// Deal with the definition 'k' of augmented to be sampled (const = FALSE).
    /// NOTE: this assume definitions are univocal and status correct.
    void helper_build_definition_order_();

    /// Set dependecies ordering for ininitializiation fo given status.
    /// At the same time assert that definition are complet and acyclic.
    void set_definition_dependency_order_for_(int status);
    /// Extract status-specific definitions row indexes.
    std::vector<size_t> get_unsorted_definition_rows_(int status, bool is_pre_process) const;
    /// Set sample map and relative usage list.
    void set_sample_map_for_(int status, size_t k);
    /// Build augmented function and inspect for constant.
    void set_function_or_identity_for_(size_t k);

    /// Deal with the definition of augmented variable with a model.
    /// NOTE: no need to initialize to generic since only used for posterior sampling.
    void helper_init_model_list_();
    /// Return dictionary of ordered models.
    hsp::str_table get_ordered_models_(int status) const;

    /// Get the definition order for current member according to its status.
    /// NOTE: if status hasn't been observed in data, build on demand.
    std::vector<size_t> helper_init_get_definition_order_(Individual &member);
    /// Assign model arguments and range. Same for function, otherwise assign proposer.
    /// NOTE: a Variable must already be allocated.
    void helper_init_attributes_(Individual &member, size_t k);
    /// Auxiliary method to assign attributes to augmented model.
    void helper_init_model_attributes_(Individual &member, size_t k);
    /// Auxiliary method to assign attributes to augmented function.
    void helper_init_func_attributes_(Individual &member, size_t k);
};
