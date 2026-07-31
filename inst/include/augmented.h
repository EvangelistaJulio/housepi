/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/nan.h"
#include "tools/random.h"
#include "tools/types.h"

#include "individual.h"
#include "variable.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  *** Augmented : struct ***  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    struct InfecProfile;
    struct DataWide;
    class Household;
    class MCMCBase;

    namespace parse
    {
        struct Entry;
    }

    /// Max resamples per member (see Augmented::restart_obs_aware).
    static constexpr size_t N_INIT_ATTEMPTS_ = 100;

    struct Augmented : public KeyModule
    {
        /// Discriminates definition rows: augmented latent vs observation levels.
        enum class Kind : std::uint8_t
        {
            Augm,
            Member,
            House,
            Status
        };

        /// Augmented latents selected for export: private ones (see str::is_hidden)
        /// and definitions that cannot be rendered are dropped.
        std::vector<size_t> list_augm_export = {};
        /// Every drawn member-level observation. No visibility filter applies: a
        /// definition name carrying a leading dot is rejected unless it is an
        /// augmented latent (see assert_valid_definition_name).
        std::vector<size_t> list_obs_member = {};
        /// House-level definition rows, in house_varia slot order. Household
        /// observations are exported as surrogate columns (see is_own_column).
        std::vector<size_t> list_obs_house = {};

        /// A surrogate column: a deterministic definition standing in for a data column.
        /// NOTE: both member- and house-level observation.
        struct Surrogate
        {
            size_t row = 0;         /// definition row in the table
            size_t slot = 0;        /// index in house_varia, household level only
            bool is_house = false;  /// reduction over members, else one value per member
            bool is_static = false; /// built from data, not drawn with the members it reads
        };
        /// Surrogate columns of both levels, ordered producer-before-consumer.
        std::vector<Surrogate> list_surrogate = {};

        /// Auxiliary block for individual augmented objects.
        struct IndivBlock
        {
            size_t idx = 0;                 /// index of status-specific definition
            std::vector<Variable> vec = {}; /// associated variables
            LogProba augm_log_lik;          /// augmented model log-likelihood
            LogProba obs_log_lik;           /// observation model log-likelihood
        };
        std::vector<IndivBlock> varia;

        /// Household-level observations, per house per observation slot.
        /// NOTE: built for both inference (data-bound) and simulation (live).
        std::vector<std::vector<Variable>> house_varia;

        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        Augmented(MCMCBase *mcmc);
        ~Augmented();

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign augmented variables to all individuals.
        void build(InfecProfile &infec_profile, const str_table &input);

        /// Set all augmented definitions to empty for all individuals.
        void clear_all();
        /// Set individual augmented definitions to empty. Used in simulations.
        void clear_definitions(Individual &m);
        /// Initialize augmented variables only if missing.
        void init_definitions_for_simulation(Individual &m);
        /// Initialize augmented dependecies and priors.
        /// NOTE: evaluation-only, never run for simulations.
        void init_for_evaluation();

        /// Request discrete augmented definition.
        void request_discrete(size_t a) { user_discrete_set_.emplace(a); }

        /// Sample observation values from their distributions.
        /// NOTE: used for simulations.
        void sample_obs(Individual &member, rnd_gen &prng);

        /// Sample the status observation, if any.
        /// NOTE: used for simulations.
        void simulate_infected_status(Household &house, Individual &member, rnd_gen &prng);

        /// Initialize surrogate columns of both levels for simulations.
        void init_house_for_simulation();
        /// Clear every member's observation Variables of a household before re-sampling.
        /// FIXME: I don't think I need this anymore
        void reset_obs_data_for_simulation(Household &house);

        /// Build the surrogate columns that stand in for data columns and point each
        /// member's read-surface at its Variable before evaluation.
        /// NOTE: only static columns are eligible, a live one has no fixed value:
        /// rejected, or skipped when exporting the data view.
        void materialize_surrogate_columns(DataWide &data_wide, bool for_export = false);

        /// Reset augmented initial values for given individual.
        void restart(Individual &member, rnd_gen &prng)
        {
            for (size_t a : list_sample(member))
                varia_of_(member).vec[a].restart(prng);
        }
        /// Draw sampled augmented from their prior for generative simulation;
        /// unlike restart it uses no proposal (simulation builds none).
        void draw_sample(Individual &member, rnd_gen &prng)
        {
            for (size_t a : list_sample(member))
                varia_of_(member).vec[a].draw_value(prng);
        }
        /// Initialize augmented initial values for all individuals.
        /// NOTE: skips automatically if no augmented to sample.
        void restart(rnd_gen &prng);
        /// Finalize all variables (see Variable::finalize).
        void finalize();

        /// Resolve individual and household ids for augmented draws and return tracked pointers.
        /// NOTE: list of draws are assumed to be unique.
        std::vector<const Individual *> init_individual_sample(
            const std::vector<std::string> &id_draws,
            const std::vector<std::string> &house_draws,
            rnd_gen &prng) const;

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Access the augmented variable of given key for a member.
        Variable &at(size_t a, const Individual &m) noexcept { return varia_of_(m).vec[a]; }
        const Variable &at(size_t a, const Individual &m) const noexcept { return varia_of_(m).vec[a]; }

        const std::vector<size_t> &list_sample(const Individual &m) const { return block_of_(m).sample; }
        const std::vector<size_t> &list_adaptive(const Individual &m) const { return block_of_(m).adaptive; }
        const std::vector<size_t> &list_draws(const Individual &m) const { return block_of_(m).draws; }
        const std::vector<size_t> &list_models(const Individual &m) const { return block_of_(m).models; }
        const std::vector<size_t> &list_obs(const Individual &m) const { return block_of_(m).obs; }

        /// Declare a read-cell for every column a definition produces, so a name
        /// standing in for no imported column resolves like any other data name.
        /// NOTE: run-scoped, an import drops the optional keys of the previous one.
        std::vector<std::string> get_produced_cols();

        /// Export column name of a surrogate column.
        const std::string &surrogate_name(const Surrogate &surr) const { return deftable_.c_name[surr.row]; }
        /// Whether a surrogate column is produced as a column of its own, rather than
        /// drawn with the observations of each infected member (see list_obs_member).
        static bool is_own_column(const Surrogate &surr) noexcept
        {
            return surr.is_house || surr.is_static;
        }
        /// How many surrogate columns are produced as columns of their own.
        size_t n_own_columns() const noexcept
        {
            size_t n = 0;
            for (const Surrogate &surr : list_surrogate)
                n += is_own_column(surr);
            return n;
        }

        /// Compute and store log-likelihood of augmented variables associated to Individual.
        /// NOTE: to be called only if any prior is defined.
        double commit_augm_log_lik(Individual &member);
        /// Compute and store log-ratio of likelihood of single augmented variable.
        double commit_augm_log_lik_ratio(Individual &member, size_t a);

        /// Whether some other variable reads 'a' as a model parameter (loc/min/max).
        /// If so, moving 'a' shifts their density and the plain value-ratio breaks.
        bool has_model_dependents(const Individual &m, size_t a) const
        {
            return block_of_(m).consumers[a].size() > 1;
        }
        /// Guarded sum of the priors that 'a' can shift (its model-consumers, incl. 'a').
        /// NOTE: out-of-support consumers contribute with -Inf.
        double sum_guarded_augm_log_lik(const Individual &m, size_t a) const;
        /// Store the log-ratio of a move of 'a' that shifts its model-consumers' densities.
        /// L_old must be sum_guarded_augm_log_lik() sampled before 'a' was proposed.
        double commit_augm_log_lik_ratio_dependent(Individual &m, size_t a, double L_old);

        /// Compute and store observation log-likelihood of observation nodes.
        /// NOTE: to be called only if any observation node is defined.
        double commit_obs_log_lik(Individual &member);
        /// Compute and store log-ratio of observation node.
        double commit_obs_log_lik_ratio(Individual &member, size_t a);

        /// Read-only individual observation log-likelihood.
        /// NOTE: an observation value is a datum, never granted valid: must be guarded.
        double eval_log_obs(const Individual &member) const
        {
            double ll = 0.0;
            for (size_t a : list_obs(member))
                ll += varia_of_(member).vec[a].guarded_log_pdf();
            return ll;
        }

        /// Throw initialization error for associated to observation model.
        [[noreturn]] void throw_obs_init_error(const Individual &member) const;

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// If there are non-constant augmented definitions.
        bool is_any_sample() const noexcept { return is_any_sample_; }
        /// If there are augmented definitions with an adaptive proposal.
        bool is_any_adaptive() const noexcept { return is_any_adaptive_; }
        /// If there are augmented definitions with a prior.
        bool is_any_model() const noexcept { return is_any_model_; }
        /// If moving given parameter may change the augmented model.
        bool param_affects_augm_model(size_t p) const noexcept { return model_param_depen_[p]; }
        /// If there are definitions for observations.
        bool is_any_obs() const noexcept { return is_any_obs_; }
        /// If a status observation is defined.
        bool is_any_status_obs() const noexcept { return !na::is(status_obs_def_); }
        /// Whether the given name is an member-level observation node.
        bool is_member_obs(const std::string &name) const noexcept;
        /// Whether the given name is a household-level observation.
        bool is_house_obs(const std::string &name) const noexcept;
        /// If moving given parameter may change the observation log-likelihood.
        bool param_affects_obs_log_lik(size_t p) const noexcept { return obs_param_depen_[p]; }
        /// If moving augmented 'a' for this member parameterizes any obs node.
        bool affects_obs(const Individual &m, size_t a) const { return !block_of_(m).obs_depen[a].empty(); }

    private:
        MCMCBase *mcmc_;               /// pointer to owning MCMCBase // UPDATE drop
        bool is_any_sample_ = false;   /// at least one augmented is to sample
        bool is_any_adaptive_ = false; /// at least one augmented has adaptive proposal
        bool is_any_model_ = false;    /// at least one augmented model is provided
        bool is_any_obs_ = false;      /// at least one observation node is provided

        size_t status_obs_def_ = na::INDEX; /// deftable row of status observation
        size_t status_obs_key_ = na::INDEX; /// key-holder slot of status observation

        /// User defined discrete parameters
        std::unordered_set<size_t> user_discrete_set_ = {};
        /// Status to definition row.
        std::unordered_map<int, std::vector<size_t>> map_status_to_def_row_ = {};
        /// Status to internal index.
        std::unordered_map<int, size_t> map_status_to_block_ = {};

        /// Parameters whose moves can change the augmented log-likelihood.
        std::vector<bool> model_param_depen_ = {};
        /// Parameters whose moves can change the observation log-likelihood.
        std::vector<bool> obs_param_depen_ = {};

        /// Auxiliary definition table.
        struct DefTable_
        {
            using col = std::vector<std::string>;
            col c_name = {}, c_status = {};
            col c_func = {}, c_fcumul = {};
            col c_value = {}, c_const = {}, c_discr = {};
            col c_prior = {}, c_min = {}, c_max = {};
            col c_proposal = {};
            std::vector<Kind> c_kind = {}; /// augm | member | house | status

            /// Clear all definitions.
            void clear();
            /// Append definitions to current table.
            void append(const str_table &input);
            /// Columns presenting definitions possibly containing sources.
            std::vector<const col *> source_cols() const
            {
                return {&c_value, &c_prior, &c_func, &c_min, &c_max};
            }
        } deftable_;

        /// Auxiliary block for status-specific priors.
        struct StatusBlock_
        {
            StatusBlock_(int status_) { status = status_; }

            int status; /// associated status

            std::vector<size_t> order = {}; /// definition order as in the deftable

            std::vector<size_t> models = {};   /// keys for augmented-model densities
            std::vector<size_t> sample = {};   /// keys for sample
            std::vector<size_t> adaptive = {}; /// keys for adaptive (subset of sample)
            std::vector<size_t> draws = {};    /// keys for draws
            /// key -> key and its transitive consumers
            std::vector<std::vector<size_t>> consumers = {};

            std::vector<size_t> obs = {};                    /// keys for observation nodes
            std::vector<std::vector<size_t>> obs_depen = {}; /// latent key -> obs nodes conditioned on it
        };
        std::vector<StatusBlock_> status_blocks_;

        /// Access IndivBlock for a member.
        IndivBlock &varia_of_(const Individual &m) noexcept { return varia[Index::of(m)]; }
        const IndivBlock &varia_of_(const Individual &m) const noexcept { return varia[Index::of(m)]; }
        /// Status-specific block governing augmented definitions of member.
        StatusBlock_ &block_at_(int status) { return status_blocks_[map_status_to_block_.at(status)]; }
        const StatusBlock_ &block_of_(const Individual &m) const { return status_blocks_[varia_of_(m).idx]; }

        /// Check augmented names and assert no reserved names are used.
        void helper_build_assert_input_(const str_table &input);

        /// Collect parameters feeding any augmented model. Infectivity profile excluded.
        /// NOTE: call once parameters are registered, i.e. from init_for_evaluation.
        /// NOTE: prior parameters may be registered after the augmented definitions.
        void set_model_param_depen_();

        /// Eventually append infectivity profile definitions.
        void helper_build_append_(InfecProfile &infec_profile);

        /// Assert augmented status/generic definitions are non degenerate.
        /// Assert status-specific definitions are only for infected status.
        /// NOTE: generic definitions are associated to MISSING status.
        void helper_build_definition_by_status_();

        /// Deal with the definition 'k' of augmented to be sampled (const = FALSE).
        /// NOTE: this assume definitions are univocal and status correct.
        void helper_build_set_allocation_order_();

        /// Build the augmented producer graph (name -> the @A names it reads across
        /// all source columns) for the given definition rows.
        str_table build_augm_graph_(const std::vector<size_t> &rows) const;

        /// Reject a read of an augmented name absent from the 'status' graph, telling
        /// apart a name declared for other statuses from a data column or a parameter
        /// of the same name (an unsuffixed name is read as augmented).
        void assert_augm_depen_declared_(const str_table &graph, int status) const;

        /// Set dependecies ordering for ininitializiation fo given status.
        /// At the same time assert that definition are complet and acyclic.
        void set_definition_dependency_order_for_(int status, int err_status);
        /// Extract status-specific definitions row indexes.
        std::vector<size_t> get_unsorted_definition_rows_(int status) const;
        /// Set sample and draws list for given status.
        void set_sample_list_(int status);
        /// Set sample and draws list for given definition.
        void set_sample_list_for_single_definition_(int status, size_t k);

        /// Inspect if constant definition.
        std::pair<bool, bool> set_sample_list_for_function_(int status, size_t k);

        /// Deal with the definition of augmented variable with a prior.
        /// NOTE: process priors only for posterior sampling, not simulations.
        void helper_init_model_list_();
        /// Return dictionary of ordered models.
        str_table get_ordered_models_(const StatusBlock_ &block) const;

        /// Deal with observation definitions.
        /// Build per-status obs lists and augmented inverse-dependency map.
        void helper_init_observation_list_();

        /// Assert status observation as correct dependency and output.
        void assert_status_obs_consistent_();

        /// Order and classify the surrogate columns of both levels into list_surrogate,
        /// and collect the household ones into list_obs_house.
        void helper_init_surrogate_list_();
        /// Whether definition row 'k' is a surrogate column materialized from data.
        bool is_surrogate_(size_t k) const noexcept;
        /// Build the Variable of a surrogate column over a house, wiring the cells it
        /// stands in for. Household level reduces over the members, member level
        /// assigns one Variable per member. Evaluation fixes the value it settles on.
        // NOTE: assume every surrogate column is registered as a data key
        void build_surrogate_column_(const Surrogate &surr, DataWide &data_wide,
                                     Household &house, bool for_evaluation);
        /// Assert a surrogate value honours what the model requires of the column.
        void assert_surrogate_datum_requirements_(const DataWide &data_wide, size_t d, double v,
                                                  const std::string &name, const Household &house) const;
        /// Whether definition row 'k' references only data and numeric sources.
        bool def_is_pure_data_(size_t k) const;
        /// Whether a parsed definition subtree references any group risk.
        static bool entry_has_group_risk_(const parse::Entry &node);

        /// Reject an augmented variable that reads a live data-bound household observation.
        void assert_no_live_house_obs_cycle_() const;
        /// Reject a definition whose prior reads an observation column.
        void assert_no_obs_entry_collision_() const;
        /// Reject a simulation override (value + prior) that reads an augmented variable
        /// or observation: it must resolve from data or parameters at draw time.
        void assert_no_sim_override_forward_dep_() const;

        /// Get the definition order for current member according to its status.
        /// NOTE: if status hasn't been observed in data, build on demand.
        std::vector<size_t> helper_materialize_definition_order_(Individual &member);
        /// Assign prior arguments and range. Same for function, otherwise assign proposer.
        /// NOTE: a Variable must already be allocated.
        /// NOTE: observations are generated by simulations.
        void helper_init_attributes_(Individual &member, size_t k, bool simulation = false);
        /// Assign attributes to augmented model.
        void helper_init_distr_attributes_(Individual &member, size_t k);
        /// Build the Variable of definition row 'k' on 'anchor'. Its level selects the
        /// operand policy: a household row reduces over the members of 'house', any other
        /// row binds every operand to 'anchor' alone and needs no household.
        void build_definition_variable_(Variable &var, size_t k,
                                        Individual &anchor, Household *house);
        /// Fold a household-level definition into a reduction over 'house': a per-member
        /// data cell unfolds into one operand each, any other operand contributes once
        /// from 'anchor'. Empty cells drop out (na.rm convention).
        /// Return whether the reduction holds any operand.
        bool bind_house_reduction_(Variable &var, const parse::Entry &root,
                                   Household &house, Individual &anchor);
    };
}
