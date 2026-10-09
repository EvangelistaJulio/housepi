/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "tools/keys.h"
#include "tools/nan.h"
#include "augmented/definition.h"
#include "data/status.h"

//////////////////////////////  augm::Scheme  //////////////////////////////

namespace housepi
{
    struct Parameter;
}

namespace housepi::augm
{
    /// Surrogate column: a deterministic definition standing in for a data column.
    ///
    /// Either a member or a household observation.
    struct Surrogate
    {
        size_t row = 0;         /// Definition row.
        size_t slot = 0;        /// Slot in `house_varia`, household level only.
        bool is_house = false;  /// Whether it reduces over members, else one value per member.
        bool is_static = false; /// Whether built from data, not drawn with the members it reads.

        /// Whether the column is produced on its own.
        /// Otherwise it is drawn with the observations of each infected member, see `list_obs_member`.
        bool is_own_column() const noexcept { return is_house || is_static; }
    };

    /// Sampling scheme of one status.
    ///
    /// - evaluation order;
    /// - what a sweep samples, adapts and draws;
    /// - densities a move recomputes.
    struct StatusScheme
    {
        std::vector<size_t> order = {};     /// Definition rows, producer before consumer, in inference.
        std::vector<size_t> sim_order = {}; /// Definition rows, producer before consumer, in simulation.

        std::vector<size_t> latent_models = {}; /// Keys carrying a latent-model density.
        std::vector<size_t> sample = {};        /// Keys to sample.
        std::vector<size_t> sim_sample = {};    /// Keys of `sample`, in simulation order.
        std::vector<size_t> adaptive = {};      /// Keys with an adaptive proposal, a subset of `sample`.
        std::vector<size_t> draws = {};         /// Keys to draw.
        /// Key -> itself and its transitive consumers, producer before consumer.
        std::vector<std::vector<size_t>> consumers = {};

        std::vector<size_t> obs_models = {}; /// Keys carrying an observation-model density, in simulation order.
        /// Key -> observation nodes transitively consuming it, itself excluded.
        std::vector<std::vector<size_t>> obs_consumers = {};

        /// Key -> whether transmission consumes it, directly or through intermediaries.
        std::vector<bool> produces_transm = {};
    };

    /// Sampling schemes of every status, and what definitions produce.
    ///
    /// Built by SchemeBuilder, fixed for a run.
    struct Scheme
    {
        std::vector<size_t> list_latent_export = {}; /// Latents to export.
        std::vector<size_t> list_obs_member = {};    /// Drawn member observations.
        std::vector<size_t> list_obs_house = {};     /// Household observation rows, in `house_varia` slot order.
        std::vector<Surrogate> list_surrogate = {};  /// Surrogate columns of both levels, producer before consumer.

        /// Status scheme at index `i`.
        const StatusScheme &at(size_t i) const noexcept { return status_schemes_[i]; }
        /// Find the index of the status scheme governing status `s`.
        /// - own scheme, if any;
        /// - else the generic one;
        /// - else na::INDEX.
        size_t index_of(StatusLvl s) const noexcept;

        /// Whether any latent is sampled.
        bool is_any_sample() const noexcept { return is_any_sample_; }
        /// Whether any latent has an adaptive proposal.
        bool is_any_adaptive() const noexcept { return is_any_adaptive_; }
        /// Whether any latent carries a model density.
        bool is_any_latent_model() const noexcept { return is_any_latent_model_; }
        /// Whether any observation carries a model density.
        bool is_any_obs_model() const noexcept { return is_any_obs_model_; }
        /// Whether any member observation column is produced, model-carrying or deterministic.
        bool is_any_obs_node() const noexcept { return !list_obs_member.empty(); }
        /// Whether definitions imply a unique household equivalence class.
        bool is_unique_eq_class() const noexcept
        {
            return is_any_sample() || is_any_latent_model() || is_any_obs_node();
        }

        /// Whether moving parameter `p` may change the latent log-likelihood.
        bool param_produces_latent_log_lik(size_t p) const noexcept { return param_produces_latent_ll_[p]; }
        /// Whether moving parameter `p` may change the observation log-likelihood.
        bool param_produces_obs_log_lik(size_t p) const noexcept { return param_produces_obs_ll_[p]; }

    private:
        friend struct SchemeBuilder;

        std::vector<StatusScheme> status_schemes_ = {}; /// Per-status schemes, then missing-observation variants.
        size_t n_base_ = 0;                             /// Number of per-status schemes, before the variants.
        std::unordered_map<StatusLvl, std::vector<size_t>> status_rows_ = {}; /// Status -> definition rows.
        std::unordered_map<StatusLvl, size_t> status_index_ = {};              /// Status -> index in `status_schemes_`.

        bool is_any_sample_ = false;       /// Whether any latent is sampled.
        bool is_any_adaptive_ = false;     /// Whether any latent has an adaptive proposal.
        bool is_any_latent_model_ = false; /// Whether any latent carries a model.
        bool is_any_obs_model_ = false;    /// Whether any observation carries a model.

        /// Parameter -> whether a latent model consumes it, so its moves shift the density.
        std::vector<bool> param_produces_latent_ll_ = {};
        /// Parameter -> whether an observation consumes it, through its model or a latent value.
        std::vector<bool> param_produces_obs_ll_ = {};
    };

    /// Build a Scheme from definitions, keyed in `keys`.
    struct SchemeBuilder
    {
        SchemeBuilder(Scheme &scheme, const DefTable &def, const KeyHolder &keys)
            : scheme_(scheme), def_(def), keys_(keys) {}

        /// Route member definition rows by status.
        /// - household observation: skipped, it needs neither key nor status.
        void register_rows();

        /// Order the definitions of every status, producer before consumer.
        ///
        /// Per status, it lists what a sweep samples, adapts and draws.
        /// Asserts that definitions are complete and acyclic.
        /// @param infected Infected statuses the data carries.
        /// @param draws_status Whether a status process draws statuses the data may not carry.
        /// @return Data names read as equivalence-class covariates.
        std::vector<std::string> order(const std::vector<StatusLvl> &infected,
                                       const Parameter &param,
                                       bool draws_status);

        /// Order and classify surrogate columns of both levels into `list_surrogate`.
        /// Household ones are also collected into `list_obs_house`.
        void init_surrogates();
        /// List observation nodes and, per key, the observations consuming it.
        void init_observations();
        /// Mark, per status, the keys that transmission names `transm` consume.
        /// - empty `transm`: no key marked.
        void init_transm(const std::vector<std::string> &transm);

        /// List latent models and, per key, the models it can shift.
        /// @note Needed for inference and evaluation only, not for simulation.
        void init_models();
        /// Drop variants appended by a previous run.
        void clear_variants() { scheme_.status_schemes_.resize(scheme_.n_base_); }
        /// Append a variant of scheme `base` leaving observations `absent` out, and return its index.
        /// A latent no sink reaches is left out too, with its draw and density.
        size_t make_variant(size_t base,
                            const std::vector<size_t> &absent,
                            const std::vector<std::string> &transm);
        /// Mark parameters feeding a latent or observation model, infectivity profile excluded.
        void init_param_producers(const KeyHolder &param);

    private:
        Scheme &scheme_;        /// Scheme being built.
        const DefTable &def_;   /// Definitions.
        const KeyHolder &keys_; /// Keys of member definitions.

        /// Sources read while ordering, and covariates collected along the way.
        struct Sources_
        {
            const Parameter &param;
            std::vector<std::string> covariates = {};
        };

        /// Status scheme of status `s`, once ordered.
        StatusScheme &at_(StatusLvl s) { return scheme_.status_schemes_[scheme_.status_index_.at(s)]; }
        /// Order the definitions of status `s`, reporting errors for `err_status`.
        void order_status_(StatusLvl s, StatusLvl err_status, Sources_ &src);
        /// List status-specific rows of `s`, completed by the generic ones.
        std::vector<size_t> unsorted_rows_(StatusLvl s) const;
        /// Add definition row `k` to the sample and draw lists of status `s`.
        void list_definition_(StatusLvl s, size_t k, Sources_ &src);
        /// Classify deterministic definition row `k` of status `s`.
        /// @return Whether to export it, and whether to draw it.
        std::pair<bool, bool> list_function_(StatusLvl s, size_t k, Sources_ &src);
    };
}
