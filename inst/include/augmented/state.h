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
        /// Index of the status scheme governing the member.
        size_t idx = 0;
        /// Latent and observation variables, by key.
        std::vector<Variable> vec = {};
        LogProba latent_log_lik;
        LogProba obs_log_lik;
    };

    /// Simulation storage of one member, bound once per run.
    struct SimMember
    {
        /// Variables per status scheme, swapped into the member's own when its status selects them.
        /// * unreachable status scheme -> no block;
        /// * key left unbound -> empty, as in State::blank;
        /// * binding error -> partial, never swapped in.
        std::vector<std::vector<Variable>> blocks = {};
        /// Binding error per status scheme, raised if reached.
        std::vector<std::exception_ptr> errors = {};
        /// Keys the member holds from its block.
        std::vector<size_t> swapped = {};
    };

    /// Mutable per-run storage of variables, moved by inference and simulation.
    struct State
    {
        /// Member variables, by individual index.
        std::vector<IndivState> varia = {};
        /// Household observations, per house and slot.
        /// * inference -> bound to data;
        /// * simulation -> drawn.
        std::vector<std::vector<Variable>> house_varia = {};

        /// Simulation storage, by individual index.
        std::vector<SimMember> sim = {};
        /// Empty variables named by key.
        std::vector<Variable> blank = {};

        /// Factors of the latent blocks, one per member with a block.
        std::vector<block_rw::Factor> block_factors = {};
        /// Factor of each individual, NA without a latent block.
        std::vector<size_t> block_factor_of = {};

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
            block_factors.clear();
            block_factor_of.clear();
        }
    };

    //////  Log-likelihood  //////

    /// Evaluate the latent log-likelihood of a member, without storing it.
    /// * value outside a support shifted by a moved parameter -> -Inf.
    inline double eval_log_latent(const StatusScheme &ss, const IndivState &st)
    {
        double ll = 0.0;
        for (size_t a : ss.latent_models)
            ll += st.vec[a].guarded_log_pdf();
        return ll;
    }
    /// Compute and store the latent log-likelihood of a member.
    inline double commit_latent_log_lik(const StatusScheme &ss, IndivState &st)
    {
        if (ss.latent_models.empty())
            return 0.;
        const double new_ll = eval_log_latent(ss, st);
        st.latent_log_lik.set(new_ll);
        return new_ll;
    }
    /// Compute and store the log-ratio of models a move of `a` touches, `a` included.
    inline double commit_latent_log_lik_ratio(const StatusScheme &ss, IndivState &st, size_t a)
    {
        double delta_ll = 0;
        // consumers[a] includes `a`: the model of `a` is counted
        for (size_t b : ss.consumers[a])
            delta_ll += st.vec[b].log_ratio();
        st.latent_log_lik.increment(delta_ll);
        return delta_ll;
    }
    /// Sum the log-densities of models `a` can shift, `a` included.
    /// * consumer out of support -> -Inf.
    inline double sum_guarded_latent_log_lik(const StatusScheme &ss, const IndivState &st, size_t a)
    {
        double ll = 0.;
        // consumers[a] includes `a`: the model of `a` is counted
        for (size_t b : ss.consumers[a])
            ll += st.vec[b].guarded_log_pdf();
        return ll;
    }
    /// Store the log-ratio of a move of `a` that shifts the densities of its model consumers.
    /// @note Take `L_old` from sum_guarded_latent_log_lik() before proposing `a`.
    inline double commit_latent_log_lik_ratio_with_consumers(const StatusScheme &ss,
                                                             IndivState &st,
                                                             size_t a,
                                                             double L_old)
    {
        // with `a` proposed, its model consumers read the moved value live
        // compare their guarded model sum against the pre-move snapshot `L_old`
        const double delta_ll = sum_guarded_latent_log_lik(ss, st, a) - L_old;
        st.latent_log_lik.increment(delta_ll);
        return delta_ll;
    }

    /// Evaluate the observation log-likelihood of a member, without storing it.
    /// * datum outside its support -> -Inf.
    inline double eval_log_obs(const StatusScheme &ss, const IndivState &st)
    {
        double ll = 0.0;
        for (size_t a : ss.obs_models)
            ll += st.vec[a].guarded_log_pdf();
        return ll;
    }
    /// Compute and store the observation log-likelihood of a member.
    inline double commit_obs_log_lik(const StatusScheme &ss, IndivState &st)
    {
        if (ss.obs_models.empty())
            return 0.;
        const double new_ll = eval_log_obs(ss, st);
        st.obs_log_lik.set(new_ll);
        return new_ll;
    }
    /// Compute and store the log-ratio of observations a move of `a` touches.
    inline double commit_obs_log_lik_ratio(const StatusScheme &ss, IndivState &st, size_t a)
    {
        // skip unless the moved latent `a` parameterizes an observation node
        if (ss.obs_consumers[a].empty())
            return 0.;

        // recompute the log-proba: an observation value is a datum
        st.obs_log_lik.set(eval_log_obs(ss, st));
        return st.obs_log_lik.ratio();
    }

    //////  Sampling  //////

    /// Initialize the sampled latents of a member for inference.
    inline void restart(const StatusScheme &ss, IndivState &st, RndGen &prng)
    {
        for (size_t a : ss.sample)
            st.vec[a].restart(prng);
    }
    /// Draw latents `keys` of a member from their models, in order.
    /// Unlike restart(), it uses no proposal.
    inline void draw_sample(const std::vector<size_t> &keys, IndivState &st, RndGen &prng)
    {
        for (size_t a : keys)
            st.vec[a].draw_value(prng);
    }
}
