/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tools/dependency.h"
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

    struct Augmented : public KeyModuleWithRequests
    {
        /// Discriminates definition rows: latent vs observation levels.
        enum class Kind : std::uint8_t
        {
            Latent,
            Member,
            House,
            Status
        };

        /// Latents selected for export.
        /// Hidden ones are not included (see str::is_hidden).
        std::vector<size_t> list_latent_export = {};
        /// Every drawn member-level observation.
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

        /// Auxiliary block for individual latent objects.
        struct IndivBlock
        {
            size_t idx = 0;                 /// index of status-specific definition
            std::vector<Variable> vec = {}; /// associated variables
            LogProba latent_log_lik;        /// latent model log-likelihood
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

        /// Declare latent keys as read by the transmission computation.
        void request_transmission(std::vector<size_t> vec);

        /// Assign latent variables to all individuals.
        void build(InfecProfile &infec_profile, const str_table &input);

        /// Set definitions to empty for every infected individual.
        /// NOTE: Used in simulations.
        void clear_all();
        /// Set individual definitions to empty.
        void clear_definitions(const Individual &m);
        /// Initialize latent variables only if missing.
        void init_definitions_for_simulation(const Individual &m);
        /// Initialize latent dependencies and models.
        /// NOTE: evaluation-only, never run for simulations.
        void init_for_evaluation();

        /// Initialize surrogate columns of both levels for simulations.
        void init_house_for_simulation();
        /// Sample the status observation, if any.
        /// UPDATE: status should not be stored in viewers.
        void simulate_infected_status(Household &house,
                                      Individual &indiv,
                                      rnd_gen &prng);
        /// Sample observation values from their distributions.
        /// NOTE: used for simulations.
        void sample_obs(const Individual &indiv, rnd_gen &prng);

        /// Build the surrogate columns that stand in for data columns and point each
        /// member's read-surface at its Variable before evaluation.
        /// NOTE: only static columns are eligible, a live one has no fixed value:
        /// rejected, or skipped when exporting the data view.
        void materialize_surrogate_columns(bool for_export = false);

        /// Reset latent initial values for given individual.
        void restart(const Individual &indiv, rnd_gen &prng)
        {
            for (size_t a : list_sample(indiv))
                varia_of_(indiv).vec[a].restart(prng);
        }
        /// Draw sampled latent from their model for generative simulation;
        /// unlike restart it uses no proposal (simulation builds none).
        void draw_sample(const Individual &indiv, rnd_gen &prng)
        {
            for (size_t a : list_sample(indiv))
                varia_of_(indiv).vec[a].draw_value(prng);
        }
        /// Initialize latent initial values for all individuals.
        /// NOTE: skips automatically if no latent to sample.
        void restart(rnd_gen &prng);
        /// Finalize all variables (see Variable::finalize).
        void finalize();

        /// Resolve individual and household ids for latent draws and return tracked pointers.
        /// NOTE: list of draws are assumed to be unique.
        std::vector<const Individual *> init_individual_sample(
            const std::vector<std::string> &id_draws,
            const std::vector<std::string> &house_draws,
            rnd_gen &prng) const;

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Access the latent variable of given key for an individual.
        Variable &at(size_t a, const Individual &m) noexcept { return varia_of_(m).vec[a]; }
        const Variable &at(size_t a, const Individual &m) const noexcept { return varia_of_(m).vec[a]; }

        const std::vector<size_t> &list_sample(const Individual &m) const { return block_of_(m).sample; }
        const std::vector<size_t> &list_adaptive(const Individual &m) const { return block_of_(m).adaptive; }
        const std::vector<size_t> &list_draws(const Individual &m) const { return block_of_(m).draws; }
        const std::vector<size_t> &list_latent_models(const Individual &m) const { return block_of_(m).latent_models; }
        const std::vector<size_t> &list_obs_models(const Individual &m) const { return block_of_(m).obs_models; }

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

        /// Number of exported columns in simulations.
        size_t n_export_cols() const noexcept;
        /// Whether each exported column has discrete-requested output type.
        /// NOTE: mixed continuous/discrete definitions are treated continuous.
        /// NOTE: it requires a simulation repetition to have run before.
        std::vector<bool> discrete_export_cols() const;

        /// Compute and store log-likelihood of latent variables associated to Individual.
        /// NOTE: to be called only if any latent model is defined.
        double commit_latent_log_lik(const Individual &indiv);
        /// Compute and store log-ratio of likelihood of single latent variable.
        double commit_latent_log_lik_ratio(const Individual &indiv, size_t a);

        /// Whether some other variable consumes 'a' as a model parameter (loc/min/max).
        /// If so, moving 'a' shifts their density and the plain value-ratio breaks.
        bool has_model_consumers(const Individual &m, size_t a) const
        {
            return block_of_(m).consumers[a].size() > 1;
        }
        /// Guarded sum of the models that 'a' can shift (its model-consumers, incl. 'a').
        /// NOTE: out-of-support consumers contribute with -Inf.
        double sum_guarded_latent_log_lik(const Individual &m, size_t a) const;
        /// Store the log-ratio of a move of 'a' that shifts its model-consumers' densities.
        /// L_old must be sum_guarded_latent_log_lik() sampled before 'a' was proposed.
        double commit_latent_log_lik_ratio_with_consumers(const Individual &m,
                                                          size_t a,
                                                          double L_old);

        /// Compute and store observation log-likelihood of observation nodes.
        /// NOTE: to be called only if any observation model is defined.
        double commit_obs_log_lik(const Individual &indiv);
        /// Compute and store log-ratio of observation node.
        double commit_obs_log_lik_ratio(const Individual &indiv, size_t a);

        /// Read-only individual observation log-likelihood.
        /// NOTE: an observation value is a datum, never granted valid: must be guarded.
        double eval_log_obs(const Individual &indiv) const
        {
            double ll = 0.0;
            for (size_t a : list_obs_models(indiv))
                ll += varia_of_(indiv).vec[a].guarded_log_pdf();
            return ll;
        }

        /// Throw initialization error for associated to observation model.
        [[noreturn]] void throw_obs_init_error(const Individual &indiv) const;

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Whether the per-individual storage covers current data.
        /// NOTE: deferred definition can leave data stale.
        bool is_allocated() const noexcept;

        /// If there are non-constant latent definitions.
        bool is_any_sample() const noexcept { return is_any_sample_; }
        /// If there are latent definitions with an adaptive proposal.
        bool is_any_adaptive() const noexcept { return is_any_adaptive_; }
        /// If any latent carries a model density.
        bool is_any_latent_model() const noexcept { return is_any_latent_model_; }
        /// If any observation carries a model density.
        bool is_any_obs_model() const noexcept { return is_any_obs_model_; }
        /// If any member-level observation column is produced, model-carrying or deterministic.
        bool is_any_obs_node() const noexcept { return !list_obs_member.empty(); }
        /// If a status observation is defined.
        bool is_any_status_obs() const noexcept { return !na::is(status_obs_def_); }

        /// Whether the given name is an member-level observation node.
        bool is_member_obs(const std::string &name) const noexcept;
        /// Whether the given name is a household-level observation.
        bool is_house_obs(const std::string &name) const noexcept;
        /// Whether definitions imply a unique household equivalence-class.
        bool is_unique_eq_class() const noexcept
        {
            return is_any_sample() || is_any_latent_model() || is_any_obs_node();
        }

        /// If moving given parameter may change the latent log-likelihood.
        bool param_produces_latent_log_lik(size_t p) const noexcept { return param_produces_latent_ll_[p]; }
        /// If moving given parameter may change the observation log-likelihood.
        bool param_produces_obs_log_lik(size_t p) const noexcept { return param_produces_obs_ll_[p]; }
        /// Whether moving latent 'a' affects transmission log-likelihood.
        bool produces_transm(const Individual &m, size_t a) const { return block_of_(m).produces_transm[a]; }
        /// If moving latent 'a' for this member parameterizes any obs node.
        bool produces_obs(const Individual &m, size_t a) const { return !block_of_(m).obs_consumers[a].empty(); }

    private:
        MCMCBase *mcmc_;                   /// pointer to owning MCMCBase // UPDATE drop
        bool is_any_sample_ = false;       /// at least one latent is to sample
        bool is_any_adaptive_ = false;     /// at least one latent has adaptive proposal
        bool is_any_latent_model_ = false; /// at least one latent carries a model
        bool is_any_obs_model_ = false;    /// at least one observation carries a model

        size_t status_obs_def_ = na::INDEX; /// deftable row of status observation
        size_t status_obs_key_ = na::INDEX; /// key-holder slot of status observation

        /// Model declared latent keys read by the transmission computation.
        std::vector<size_t> user_transm_vec_ = {};
        /// Status to definition row.
        std::unordered_map<StatusLevel, std::vector<size_t>> map_status_to_def_row_ = {};
        /// Status to internal index.
        std::unordered_map<StatusLevel, size_t> map_status_to_block_ = {};

        /// Parameter -> whether a latent model consumes it, so its moves shift the density.
        std::vector<bool> param_produces_latent_ll_ = {};
        /// Parameter -> whether an observation consumes it, through its model or a latent value.
        std::vector<bool> param_produces_obs_ll_ = {};

        /// Auxiliary definition table.
        struct DefTable_
        {
            using col = std::vector<std::string>;
            col c_name = {}, c_status = {};
            col c_func = {}, c_fcumul = {};
            col c_value = {}, c_const = {}, c_discr = {};
            col c_model = {}, c_min = {}, c_max = {};
            col c_proposal = {};
            /// augm | member | house | status
            std::vector<Kind> c_kind = {};

            /// Whether the table is empty.
            bool empty() const noexcept { return c_name.empty(); }
            /// Clear all definitions.
            void clear();
            /// Append definitions to current table.
            void append(const str_table &input);
            /// Columns presenting definitions possibly containing sources.
            std::vector<const col *> source_cols() const
            {
                return {&c_value, &c_model, &c_func, &c_min, &c_max};
            }
        } deftable_;

        /// Auxiliary block for status-specific models.
        struct StatusBlock_
        {
            std::vector<size_t> order = {}; /// definition order as in the deftable

            std::vector<size_t> latent_models = {}; /// keys for latent-model densities
            std::vector<size_t> sample = {};        /// keys for sample
            std::vector<size_t> adaptive = {};      /// keys for adaptive (subset of sample)
            std::vector<size_t> draws = {};         /// keys for draws
            /// key -> the key and its transitive consumers, producer before consumer
            std::vector<std::vector<size_t>> consumers = {};

            std::vector<size_t> obs_models = {}; /// keys for observation-model densities
            /// key -> the obs nodes that transitively consume it, itself excluded
            std::vector<std::vector<size_t>> obs_consumers = {};

            /// key -> whether transmission consumes it, directly or through intermediaries
            std::vector<bool> produces_transm = {};
        };
        std::vector<StatusBlock_> status_blocks_;

        /// Access IndivBlock for an individual.
        IndivBlock &varia_of_(const Individual &m) noexcept { return varia[Index::of(m)]; }
        const IndivBlock &varia_of_(const Individual &m) const noexcept { return varia[Index::of(m)]; }
        /// Status-specific block governing latent definitions of member.
        StatusBlock_ &block_at_(StatusLevel s) { return status_blocks_[map_status_to_block_.at(s)]; }
        const StatusBlock_ &block_of_(const Individual &m) const { return status_blocks_[varia_of_(m).idx]; }

        // Restore internal list and flags definitions.
        void build_restore_internal_();
        /// Check latent names and assert no reserved names are used.
        void build_assert_input_(const str_table &input);
        /// Eventually append infectivity profile definitions.
        void build_append_profile_(InfecProfile &infec_profile);

        /// Collect parameters feeding any latent model. Infectivity profile excluded.
        /// NOTE: call once parameters are registered, i.e. from init_for_evaluation.
        /// NOTE: model parameters may be registered after the latent definitions.
        void set_param_producer_lists_();

        /// Assert latent status/generic definitions are non degenerate.
        /// Assert status-specific definitions are only for infected status.
        /// NOTE: generic definitions are associated to MISSING status.
        /// NOTE: assertion of univocal name-status is done in R
        void build_definition_by_status_();

        /// Deal with the definition 'k' of latent to be sampled (const = FALSE).
        /// NOTE: this assume definitions are univocal and status correct.
        void build_set_allocation_order_();

        /// Build the latent producer graph (name -> the @A names it reads across
        /// all source columns) for the given definition rows.
        dpdcy::prod_graph build_producer_graph_(const std::vector<size_t> &rows) const;

        /// Reject a read of a latent name absent from the 'status' graph.
        /// Tell apart a name declared for other statuses,
        /// from an unsuffixed name is read as latent.
        void assert_producers_declared_(const dpdcy::prod_graph &graph, StatusLevel s) const;

        /// Order the definitions of given status, producer before consumer.
        /// At the same time assert that definitions are complete and acyclic.
        void set_definition_order_for_(StatusLevel s, StatusLevel err_status);

        /// Extract status-specific definitions row indexes.
        std::vector<size_t> get_unsorted_definition_rows_(StatusLevel s) const;
        /// Set sample and draws list for given status.
        void set_sample_list_(StatusLevel s);
        /// Set sample and draws list for given definition.
        void set_sample_list_for_single_definition_(StatusLevel s, size_t k);

        /// Inspect if constant definition.
        std::pair<bool, bool> set_sample_list_for_function_(StatusLevel s, size_t k);

        /// Deal with the definition of latent variable with a model.
        /// Fills StatusBlock_::consumers, i.e. keyed by producer.
        /// NOTE: process models only for posterior sampling, not simulations.
        void init_model_list_();
        /// Model name -> its transitive consumers among the other model nodes.
        str_table get_model_consumers_(const StatusBlock_ &block) const;

        /// Deal with observation definitions.
        /// Fills StatusBlock_::obs_consumers, i.e. keyed by producer.
        void init_observation_list_();

        /// Mark the latent keys the transmission side consumes, per status.
        /// NOTE: follows model request.
        void init_transm_producers_();

        /// Assert status observation as correct dependency and output.
        void assert_status_obs_consistent_();

        /// Order and classify the surrogate columns of both levels into list_surrogate,
        /// and collect the household ones into list_obs_house.
        void init_surrogate_list_();
        /// Build the Variable of a surrogate column over a house, wiring the cells it
        /// stands in for. Household level reduces over the members, member level
        /// assigns one Variable per member. Evaluation fixes the value it settles on.
        // NOTE: assume every surrogate column is registered as a data key
        void build_surrogate_column_(const Surrogate &surr,
                                     const Household &house,
                                     bool for_evaluation);
        /// Assert a surrogate datum matches the model requests for the column it fills.
        void assert_surrogate_datum_requests_(size_t d,
                                              const Variable &var,
                                              const std::string &context) const;

        /// Whether member variable filling one exported column is discrete.
        bool is_member_col_discrete_(size_t a) const;
        /// Whether household variable filling one exported column is discrete.
        bool is_house_col_discrete_(size_t slot) const;
        /// Whether any status-specific declaration is declared continuous.
        bool has_continuous_declaration_(const std::string &name) const;

        /// Whether definition row 'k' references only data and numeric sources.
        bool def_is_pure_data_(size_t k) const;
        /// Whether a parsed definition subtree references any group risk.
        static bool entry_has_group_risk_(const parse::Entry &node);

        /// Reject a latent variable that reads a live data-bound household observation.
        void assert_no_live_house_obs_cycle_() const;
        /// Reject a definition whose model reads an observation column.
        void assert_no_obs_entry_collision_() const;
        /// Reject a simulation override (value + model) that reads a latent variable
        /// or observation: it must resolve from data or parameters at draw time.
        void assert_no_sim_override_latent_read_() const;

        /// Get the definition order for current member according to its status.
        /// NOTE: if status hasn't been observed in data, build on demand.
        std::vector<size_t> helper_materialize_definition_order_(const Individual &indiv);
        /// Assign model arguments and domain. Same for function, otherwise assign proposer.
        /// NOTE: a Variable must already be allocated.
        /// NOTE: observations are generated by simulations.
        void helper_init_attributes_(const Individual &indiv,
                                     size_t k,
                                     bool simulation = false);
        /// Assign attributes to latent model.
        void helper_init_distr_attributes_(const Individual &indiv, size_t k);
        /// Build the Variable of definition row 'k' on 'anchor'. Its level selects the
        /// operand policy: a household row reduces over the members of 'house', any other
        /// row binds every operand to 'anchor' alone and needs no household.
        void build_definition_variable_(Variable &var,
                                        size_t k,
                                        const Individual &anchor,
                                        const Household *house);
        /// Fold a household-level definition into a reduction over 'house': a per-member
        /// data cell unfolds into one operand each, any other operand contributes once
        /// from 'anchor'. Empty cells drop out (na.rm convention).
        /// Return whether the reduction holds any operand.
        bool bind_house_reduction_(Variable &var,
                                   const parse::Entry &root,
                                   const Household &house,
                                   const Individual &anchor);
    };
}
