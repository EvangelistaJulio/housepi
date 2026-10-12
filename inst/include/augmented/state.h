/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <exception>
#include <vector>

#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/random.h"
#include "augmented/scheme.h"
#include "data/individual.h"
#include "variable.h"

//////////////////////////////  augm::State  //////////////////////////////

namespace housepi::augm
{
    /// Latent and observation variables of one member, with their log-likelihoods.
    struct IndivState
    {
        size_t idx = 0;                 /// Index of the status scheme governing the member.
        std::vector<Variable> vec = {}; /// Latent and observation variables, by key.
        LogProba latent_log_lik;        /// Log-likelihood of latent models.
        LogProba obs_log_lik;           /// Log-likelihood of observation models.
    };

    /// Simulation storage of one member, bound once per run.
    struct SimMember
    {
        /// Variables per status scheme, swapped into the member's own when its status selects them.
        /// - unreachable status scheme: no block;
        /// - key left unbound: empty, as in State::blank;
        /// - binding error: partial, never swapped in.
        std::vector<std::vector<Variable>> blocks = {};
        std::vector<std::exception_ptr> errors = {}; /// Binding error per status scheme, raised if reached.
        std::vector<size_t> swapped = {};            /// Keys the member holds from its block.
    };

    /// Mutable per-run storage of variables, moved by inference and simulation.
    struct State
    {
        std::vector<IndivState> varia = {}; /// Member variables, by individual index.
        /// Household observations, per house and slot.
        /// - inference: bound to data;
        /// - simulation: live.
        std::vector<std::vector<Variable>> house_varia = {};

        std::vector<SimMember> sim = {};  /// Simulation storage, by individual index.
        std::vector<Variable> blank = {}; /// Empty variables named by key.

        /// State of member `m`.
        IndivState &of(const Individual &m) noexcept { return varia[Index::of(m)]; }
        const IndivState &of(const Individual &m) const noexcept { return varia[Index::of(m)]; }

        /// Allocate variables of `n_indiv` members for every key, named after their level.
        void allocate(size_t n_indiv,
                      const KeyHolder &keys,
                      const DefTable &def);
        /// Drop every variable.
        void clear() noexcept
        {
            varia.clear();
            house_varia.clear();
            sim.clear();
            blank.clear();
        }
    };

    //////  Log-likelihood  //////

    /// Compute and store the latent-model log-likelihood of a member.
    /// - value outside a support shifted by a moved parameter: -Inf.
    double commit_latent_log_lik(const StatusScheme &ss,
                                 IndivState &st);
    /// Compute and store the log-ratio of models a move of `a` touches, `a` included.
    double commit_latent_log_lik_ratio(const StatusScheme &ss,
                                       IndivState &st,
                                       size_t a);
    /// Sum the log-densities of models `a` can shift, `a` included.
    /// - consumer out of support: -Inf.
    double sum_guarded_latent_log_lik(const StatusScheme &ss,
                                      const IndivState &st,
                                      size_t a);
    /// Store the log-ratio of a move of `a` that shifts the densities of its model consumers.
    /// @note Take `L_old` from sum_guarded_latent_log_lik() before proposing `a`.
    double commit_latent_log_lik_ratio_with_consumers(const StatusScheme &ss,
                                                      IndivState &st,
                                                      size_t a,
                                                      double L_old);

    /// Evaluate the observation log-likelihood of a member, without storing it.
    /// - datum outside its support: -Inf.
    double eval_log_obs(const StatusScheme &ss,
                        const IndivState &st);
    /// Compute and store the observation log-likelihood of a member.
    double commit_obs_log_lik(const StatusScheme &ss,
                              IndivState &st);
    /// Compute and store the log-ratio of observations a move of `a` touches.
    double commit_obs_log_lik_ratio(const StatusScheme &ss,
                                    IndivState &st, size_t a);

    //////  Sampling  //////

    /// Initialize the sampled latents of a member for inference.
    void restart(const StatusScheme &ss, IndivState &st, rnd_gen &prng);
    /// Draw latents `keys` of a member from their models, in order.
    /// Unlike restart(), it uses no proposal.
    void draw_sample(const std::vector<size_t> &keys, IndivState &st, rnd_gen &prng);
}
