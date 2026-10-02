/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
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
    struct InfecProfile;
    namespace data
    {
        struct Wide;
    }
    class Household;
    class MCMCBase;

    /// Latent variables and observations of a run.
    /// Based on @cite tanner_calculation_1987
    ///
    /// - definitions, stored;
    /// - sampling schemes, built from definitions;
    /// - member and household variables, bound to definitions.
    struct Augmented : public KeyModuleWithRequests
    {
        augm::State state; /// Variables moved by inference and simulation.

        //////  Constructor  //////

        Augmented(MCMCBase *mcmc);
        ~Augmented();

        //////  Setter  //////

        /// Declare latent keys `vec` as read by transmission.
        void request_transmission(std::vector<size_t> vec);

        /// Build definitions from `input` and the infectivity profile, then their sampling schemes.
        void build(InfecProfile &infec_profile, const str_table &input);

        /// Empty the definitions of every infected individual, for simulation.
        void clear_all();
        /// Empty the definitions of member `m`.
        void clear_definitions(const Individual &m);
        /// Initialize the missing latents of member `m`, for simulation.
        void init_definitions_for_simulation(const Individual &m);
        /// Initialize latent dependencies and models.
        /// @note Evaluation only, never for simulation.
        void init_for_evaluation();

        /// Initialize surrogate columns of both levels, for simulation.
        void init_house_for_simulation();
        /// Draw the status observation of `indiv`, if any.
        /// @todo Status should not be stored in viewers.
        void simulate_infected_status(Household &house,
                                      Individual &indiv,
                                      rnd_gen &prng);
        /// Draw the observations of `indiv` from their models, for simulation.
        void sample_obs(const Individual &indiv, rnd_gen &prng)
        {
            augm::sample_obs(scheme_of_(indiv), state.of(indiv), prng);
        }

        /// Build surrogate columns standing in for data columns, and point member read-surfaces at them.
        ///
        /// Only static columns are eligible: a live one has no fixed value.
        /// - evaluation: live column rejected;
        /// - `for_export`: live column skipped.
        void materialize_surrogate_columns(bool for_export = false);

        /// Initialize the latents of `indiv`, for inference.
        void restart(const Individual &indiv, rnd_gen &prng)
        {
            augm::restart(scheme_of_(indiv), state.of(indiv), prng);
        }
        /// Initialize the latents of every individual, for inference.
        /// - no latent to sample: skipped.
        void restart(rnd_gen &prng);

        /// Draw the sampled latents of `indiv` from their models, for simulation.
        /// Unlike restart(), it uses no proposal.
        void draw_sample(const Individual &indiv, rnd_gen &prng)
        {
            augm::draw_sample(scheme_of_(indiv), state.of(indiv), prng);
        }
        /// Draw the sampled latents of every individual from their models, for simulation.
        /// Unlike restart(), it uses no proposal.
        /// - no latent to sample: skipped.
        void draw_sample(rnd_gen &prng);

        /// Finalize every variable, see Variable::finalize().
        void finalize();

        //////  Getter  //////

        /// Latent and observation definitions.
        const augm::DefTable &def() const noexcept { return def_; }
        /// Sampling schemes of every status, and what definitions produce.
        const augm::Scheme &scheme() const noexcept { return scheme_; }

        /// Variable of latent `a` for member `m`.
        Variable &at(size_t a, const Individual &m) noexcept { return state.of(m).vec[a]; }
        const Variable &at(size_t a, const Individual &m) const noexcept { return state.of(m).vec[a]; }

        /// Keys to sample for member `m`.
        const std::vector<size_t> &list_sample(const Individual &m) const { return scheme_of_(m).sample; }
        /// Keys with an adaptive proposal for member `m`.
        const std::vector<size_t> &list_adaptive(const Individual &m) const { return scheme_of_(m).adaptive; }
        /// Keys to draw for member `m`.
        const std::vector<size_t> &list_draws(const Individual &m) const { return scheme_of_(m).draws; }
        /// Keys carrying a latent-model density for member `m`.
        const std::vector<size_t> &list_latent_models(const Individual &m) const { return scheme_of_(m).latent_models; }
        /// Keys carrying an observation-model density for member `m`.
        /// - missing datum: marginalized out, see init_live_schemes_().
        const std::vector<size_t> &list_obs_models(const Individual &m) const { return scheme_of_(m).obs_models; }

        /// List the columns definitions produce: drawn member observations and surrogate columns.
        std::vector<std::string> get_produced_cols();

        /// List the columns drawn, or computed from a draw, in simulation.
        std::vector<std::string> get_live_cols() const;

        /// Number of exported columns in simulation.
        size_t n_export_cols() const noexcept;
        /// Whether each exported column requests a discrete output type.
        /// - mixed continuous and discrete definitions: continuous.
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
        /// - deferred definitions: storage may be stale.
        bool is_allocated() const noexcept;

        /// Whether moving latent `a` of member `m` changes the transmission log-likelihood.
        bool produces_transm(const Individual &m, size_t a) const { return scheme_of_(m).produces_transm[a]; }
        /// Whether moving latent `a` of member `m` changes an observation node.
        bool produces_obs(const Individual &m, size_t a) const { return !scheme_of_(m).obs_consumers[a].empty(); }

    private:
        MCMCBase *mcmc_; /// Owning MCMCBase. @todo Drop.

        std::vector<size_t> user_transm_vec_ = {}; /// Model-declared latent keys read by transmission.

        augm::DefTable def_;  /// Latent and observation definitions.
        augm::Scheme scheme_; /// Sampling schemes built from definitions.

        /// Status scheme governing the definitions of member `m`.
        const augm::StatusScheme &scheme_of_(const Individual &m) const { return scheme_.at(state.of(m).idx); }
        /// Build a Binder of definitions over current data.
        augm::Binder binder_();

        /// Reset the internal lists and flags of definitions.
        void build_restore_internal_();
        /// Append infectivity-profile definitions, if any.
        void build_append_profile_(InfecProfile &infec_profile);

        /// Point every individual at a status scheme holding what its run reads.
        /// - observations whose datum is present;
        /// - definitions still reaching a sink.
        ///
        /// A datum absent at build stays absent for the run: the likelihood never reads it.
        /// @note Evaluation only: a simulation reads every definition it produces.
        void init_live_schemes_(augm::SchemeBuilder &builder);

        /// List names the transmission reads: model-declared keys and the latent infectivity profile.
        std::vector<std::string> transm_names_() const;

        /// Whether the member variable filling exported column `a` is discrete.
        bool is_member_col_discrete_(size_t a) const;
        /// Whether the household variable filling exported slot `slot` is discrete.
        bool is_house_col_discrete_(size_t slot) const;

        /// Reject a member observation filling an every-member column only for infected members.
        void assert_every_member_cols_() const;

        /// Point member `indiv` at the scheme of its status, and return its definition order.
        /// - status the data does not carry: generic scheme.
        const std::vector<size_t> &helper_materialize_definition_order_(const Individual &m);
    };
}
