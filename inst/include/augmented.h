/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "tools/random.h"
#include "augmented/binder.h"
#include "augmented/definition.h"
#include "augmented/scheme.h"
#include "augmented/state.h"
#include "data/individual.h"
#include "variable.h"

//////////////////////////////  Augmented  //////////////////////////////

namespace housepi
{
    struct Data;
    struct InfecProfile;
    struct Parameter;
    struct StatusProcess;
    struct Translator;
    class Household;

    /// Latent variables and observations of a run.
    /// Based on @cite tanner_calculation_1987
    ///
    /// * definitions, stored;
    /// * sampling schemes, built from definitions;
    /// * member and household variables, bound to definitions.
    struct Augmented : public KeyHolderWithRequests
    {
        /// Variables moved by inference and simulation.
        augm::State state;

        //////  Constructor  //////

        /// Register model-declared latent `keys`, and bind the modules of its chain.
        /// @note The Translator may be bound before it is constructed: it is first read at build.
        Augmented(std::vector<std::string> keys,
                  Data &data,
                  const Parameter &param,
                  InfecProfile &infec_profile,
                  const StatusProcess &status_process,
                  const Translator &translator);
        ~Augmented();

        //////  Setter  //////

        /// Declare latent keys `vec` as read by transmission.
        void request_transmission(std::vector<size_t> vec);

        /// Build definitions from `input` and the infectivity profile, then their sampling schemes.
        /// * status process defined -> every status it draws needs a scheme.
        void build(const str_table &input);

        /// Empty the variables of every infected member, between simulated repetitions.
        /// * variable from a block -> returned to it;
        /// * value transmission set -> cleared;
        /// * static deterministic observation -> kept.
        /// * no definitions allocated -> no-op.
        void clear_all();
        /// Hand member `m` the variables of its status scheme for simulation.
        /// * variable transmission set -> kept;
        /// * binding error of the scheme, kept by init_members_for_simulation() -> raised.
        void init_definitions_for_simulation(const Individual &m);
        /// Initialize latent dependencies and models.
        /// @note Evaluation only, never for simulation.
        void init_for_evaluation();

        /// Initialize deterministic observations of both levels, for simulation.
        void init_house_for_simulation();
        /// Bind the member variables a simulation reads, once per run.
        ///
        /// * one block per status scheme the member can reach, handed out by init_definitions_for_simulation();
        /// * reachable status: the observed one, and every one the status process can draw;
        /// * each block finalized while bound -> constants folded per status scheme.
        ///
        /// A binding error is kept, and raised only if the member reaches what failed.
        /// Distinct warnings are raised once each.
        /// @note Call after init_house_for_simulation(), and release_simulation() when the run ends.
        void init_members_for_simulation();
        /// Empty every member variable and drop the blocks, at the end of a simulation run.
        void release_simulation();
        /// Draw the observations of `indiv` from their models, for simulation.
        /// * each draw lands in its Variable, where data readers read it;
        /// * the data column, if any, is never overwritten.
        void sample_obs(const Individual &indiv, RndGen &prng)
        {
            augm::draw_sample(scheme_of_(indiv).obs_models, state.of(indiv), prng);
        }

        /// Build columns of deterministic observations, and point member read-surfaces at them.
        ///
        /// Only static columns are eligible: a drawn one has no fixed value.
        /// * evaluation -> drawn column rejected;
        /// * `for_export` -> drawn column skipped.
        void materialize_deterministic_obs(bool for_export = false);

        /// Initialize the latents of `indiv`, for inference.
        void restart(const Individual &indiv, RndGen &prng)
        {
            augm::restart(scheme_of_(indiv), state.of(indiv), prng);
        }
        /// Initialize the latents of every individual, for inference.
        /// * no latent to sample -> skipped.
        void restart(RndGen &prng);

        /// Draw the sampled latents of `indiv` from their models, for simulation.
        /// Unlike restart(), it uses no proposal.
        void draw_sample_for_simulation(const Individual &indiv, RndGen &prng)
        {
            augm::draw_sample(scheme_of_(indiv).sim_sample, state.of(indiv), prng);
        }
        /// Draw the sampled latents of every individual from their models, for evaluation.
        /// Unlike restart(), it uses no proposal.
        /// * no latent to sample -> skipped;
        /// * latent pinned by pin_values() -> set to its value, read by the latents after it.
        void draw_sample(RndGen &prng);
        /// Pin latents to the columns of `table`, for draw_sample().
        /// * column -> latent of the same name, one value per member in dataset order;
        /// * `NA` cell, latent without a column -> not pinned.
        /// @note Pins set values only: sampling schemes and equivalence classes stay as built.
        void pin_values(const dbl_table &table);

        /// Finalize every variable, see Variable::finalize().
        void finalize();

        //////  Getter  //////

        /// Latent and observation definitions.
        const augm::DefTable &def() const noexcept { return def_; }
        /// Sampling schemes of every status, and what definitions produce.
        const augm::Scheme &scheme() const noexcept { return scheme_; }
        /// Latents the model generates through transmission.
        std::vector<std::string> transmission_keys() const { return consult(user_transm_vec_); }

        /// Variable of latent `a` for member `m`.
        Variable &at(size_t a, const Individual &m) noexcept { return state.of(m).vec[a]; }
        const Variable &at(size_t a, const Individual &m) const noexcept { return state.of(m).vec[a]; }
        /// Factor of the latent block of member `m`, nullptr without one.
        block_rw::Factor *block_factor(const Individual &m) noexcept
        {
            const size_t b = state.block_factor_of[Index::of(m)];
            return na::is(b) ? nullptr : &state.block_factors[b];
        }
        /// Keys to sample for member `m`.
        const std::vector<size_t> &list_sample(const Individual &m) const { return scheme_of_(m).sample; }
        /// Keys with an adaptive proposal for member `m`.
        const std::vector<size_t> &list_adaptive(const Individual &m) const { return scheme_of_(m).adaptive; }
        /// Keys of the latent block of member `m`.
        const std::vector<size_t> &list_block(const Individual &m) const { return scheme_of_(m).block; }
        /// Keys to sample for member `m`, outside its latent block.
        const std::vector<size_t> &list_sample_alone(const Individual &m) const { return scheme_of_(m).sample_alone; }
        /// Terms a move of the latent block of member `m` touches.
        const LogLikTerms &block_terms(const Individual &m) const { return scheme_of_(m).block_terms; }
        /// Keys to draw for member `m`.
        const std::vector<size_t> &list_draws(const Individual &m) const { return scheme_of_(m).draws; }
        /// Keys carrying a latent-model density for member `m`.
        const std::vector<size_t> &list_latent_models(const Individual &m) const { return scheme_of_(m).latent_models; }
        /// Keys carrying an observation-model density for member `m`.
        /// * missing datum -> marginalized out, see init_live_schemes_().
        const std::vector<size_t> &list_obs_models(const Individual &m) const { return scheme_of_(m).obs_models; }

        /// List the columns definitions produce: drawn member observations and deterministic observations.
        std::vector<std::string> get_produced_cols();

        /// List the columns drawn, or computed from a draw, in simulation.
        std::vector<std::string> get_drawn_cols() const;

        /// Number of exported columns in simulation.
        size_t n_export_cols() const noexcept;
        /// Whether each exported column requests a discrete output type.
        /// * mixed continuous and discrete definitions -> continuous.
        /// @note Run a simulation repetition first.
        std::vector<bool> discrete_export_cols() const;

        /// Compute and store the latent log-likelihood of member `m`.
        /// @note Call only if a latent model is defined.
        double commit_latent_log_lik(const Individual &m)
        {
            return augm::commit_latent_log_lik(scheme_of_(m), state.of(m));
        }
        /// Compute and store the latent log-ratio of a move of `a` for member `m`.
        double commit_latent_log_lik_ratio(const Individual &m, size_t a)
        {
            return augm::commit_latent_log_lik_ratio(scheme_of_(m), state.of(m), a);
        }

        /// Whether another variable reads `a` as a model argument (loc, min, max).
        /// If so, moving `a` shifts their density and the plain value ratio breaks.
        bool has_model_consumers(const Individual &m, size_t a) const
        {
            return scheme_of_(m).consumers[a].size() > 1;
        }
        /// Sum the log-densities of models `a` can shift, see augm::sum_guarded_latent_log_lik().
        double sum_guarded_latent_log_lik(const Individual &m, size_t a) const
        {
            return augm::sum_guarded_latent_log_lik(scheme_of_(m), state.of(m), a);
        }
        /// Store the log-ratio of a move of `a` that shifts the densities of its model consumers.
        /// @note Take `L_old` from sum_guarded_latent_log_lik() before proposing `a`.
        double commit_latent_log_lik_ratio_with_consumers(const Individual &m, size_t a, double L_old)
        {
            return augm::commit_latent_log_lik_ratio_with_consumers(scheme_of_(m), state.of(m), a, L_old);
        }

        /// Compute and store the observation log-likelihood of member `m`.
        /// @note Call only if an observation model is defined.
        double commit_obs_log_lik(const Individual &m)
        {
            return augm::commit_obs_log_lik(scheme_of_(m), state.of(m));
        }
        /// Compute and store the observation log-ratio of a move of `a` for member `m`.
        double commit_obs_log_lik_ratio(const Individual &m, size_t a)
        {
            return augm::commit_obs_log_lik_ratio(scheme_of_(m), state.of(m), a);
        }
        /// Evaluate the observation log-likelihood of member `m`, without storing it.
        double eval_log_obs(const Individual &m) const
        {
            return augm::eval_log_obs(scheme_of_(m), state.of(m));
        }

        /// Throw the initialization error caused by the observations of `indiv`.
        [[noreturn]] void throw_obs_init_error(const Individual &m) const;

        //////  Flags  //////

        /// Whether per-individual storage covers current data.
        /// * deferred definitions -> storage may be stale.
        bool is_allocated() const noexcept;

        /// Whether moving latent `a` of member `m` changes the transmission log-likelihood.
        bool produces_transm(const Individual &m, size_t a) const { return scheme_of_(m).produces_transm[a]; }
        /// Whether moving latent `a` of member `m` changes an observation node.
        bool produces_obs(const Individual &m, size_t a) const { return !scheme_of_(m).obs_consumers[a].empty(); }

    private:
        /// Data tables, read-cells and equivalence-class keys.
        /// * build -> writes equivalence-class keys;
        /// * deterministic observations, simulation init -> wire read-cells.
        Data &data_;
        /// Parameters read by definitions.
        const Parameter &param_;
        /// Infectivity profile.
        /// * latent -> definitions appended, latent key set at build.
        InfecProfile &infec_profile_;
        /// Status process, whether defined at build, statuses a member can draw in simulation.
        const StatusProcess &status_process_;
        /// Translator of definition sources.
        const Translator &translator_;

        /// Model-declared latent keys generated and read by transmission.
        std::vector<size_t> user_transm_vec_ = {};
        /// Key -> whether simulation binds it for the whole run.
        std::vector<bool> sim_fixed_ = {};
        /// Key -> pinned value per member, empty when not pinned.
        std::vector<std::vector<double>> pins_ = {};

        /// Latent and observation definitions.
        augm::DefTable def_;
        /// Sampling schemes built from definitions.
        augm::Scheme scheme_;

        /// Status scheme governing the definitions of member `m`.
        const augm::StatusScheme &scheme_of_(const Individual &m) const { return scheme_.at(state.of(m).idx); }
        /// Build a Binder of definitions over current data.
        augm::Binder binder_();

        /// Reset the internal lists and flags of definitions.
        void build_restore_internal_();
        /// Append infectivity-profile definitions, if any.
        void build_append_profile_();

        /// Point every individual at a status scheme holding what its run reads.
        /// * observations whose datum is present;
        /// * definitions still reaching a sink.
        ///
        /// A datum absent at build stays absent for the run: the likelihood never reads it.
        /// @note Evaluation only: a simulation reads every definition it produces.
        void init_live_schemes_(augm::SchemeBuilder &builder);
        /// Set the latent block of every status scheme an infected member follows, and its members' factors.
        void init_latent_blocks_(augm::SchemeBuilder &builder);
        /// List the sampled latents of member `indiv` moving together, producer before consumer.
        /// * proposal declared, or discrete -> left out, it moves alone;
        /// * bounds or location shift moving with another block latent, through functions too -> moves alone.
        std::vector<size_t> list_block_keys_(const Individual &indiv) const;

        /// List names the transmission reads: model-declared keys and the latent infectivity profile.
        std::vector<std::string> transm_names_() const;

        /// Whether the variables filling an exported column are discrete, at least one carrying a value.
        /// * member column `i` -> a key, never discrete when declared continuous;
        /// * household column `i` -> a slot.
        bool is_col_discrete_(size_t i, bool is_house) const;

        /// Reject a member observation filling an every-member column only for infected members.
        void assert_every_member_cols_() const;

        /// Bind the blocks of member `m`, one per status scheme it can reach.
        void bind_simulation_blocks_(const augm::Binder &binder, const Individual &m);

        /// Point member `m` at the scheme of its status, and return its definition order for the run.
        /// * status the data does not carry -> generic scheme.
        const std::vector<size_t> &helper_materialize_definition_order_(const Individual &m, bool simulation);
    };
}
