/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/logproba.h"
#include "tools/nan.h"
#include "augmented.h"
#include "data.h"
#include "infec_profile.h"
#include "mcmc/posterior.h"
#include "mcmc/settings.h"
#include "parameter.h"
#include "proposal.h"
#include "rates.h"
#include "status_process.h"
#include "translator.h"

//////////////////////////////  Workspace  //////////////////////////////

namespace housepi
{
    class ModelBase;

    /// Log-likelihood terms a parameter move recomputes.
    struct ParamRoute
    {
        /// Terms the move shifts.
        LogLikTerms terms = {};
        /// Equivalence classes reading the parameter.
        /// * read by more than rates, or by every class -> empty, every class recomputed.
        std::vector<size_t> classes = {};
    };

    /// Parameters moving together, one block per baseline and the rates coupled to it.
    struct ParamBlock
    {
        /// Parameters, baseline first.
        std::vector<size_t> keys = {};
        /// Union of the routes of the members.
        ParamRoute route = {};
        /// Jump factor and running moments of the coordinates.
        block_rw::Factor factor = {};
    };

    /// Definitions of a run, parsed from the imported tables, read by every chain.
    /// * built once per import, never written while sampling;
    /// * a workspace builds its modules from them, slot by slot, see Workspace::build().
    ///
    /// @note Replaying every slot in `Slot` order rebuilds a workspace only after the wall push:
    /// before it, the latent slot may wait on incomplete definitions.
    /// @note A rebuilt workspace loses members set outside these slots, unless the rebuild sets them.
    struct Definitions
    {
        /// Slots a workspace builds from, in replay order.
        enum class Slot : std::uint8_t
        {
            Param,
            GroupRate,
            ContactRate,
            InfecProfile,
            Latent,
            ParamPins,
            LatentPins
        };

        /// Parameter definitions.
        param::DefTable param;
        /// Group-rate rules, resolved per member at build.
        str_table group_rate;
        /// Contact-rate rules, resolved per household at build.
        str_table contact_rate;
        /// Infectivity-profile definitions.
        str_table infec_profile;
        /// Latent and observation definitions, then the latent infectivity-profile rows.
        augm::DefTable latent;
        /// Status drawn per member at infection, in simulation.
        StatusProcess status_process;
        /// Pinned parameter values, one per simulation repetition, by name.
        dbl_table param_pins;
        /// Pinned latent values, one per member, by name.
        dbl_table latent_pins;

        /// Set latent definitions to `input`, then the rows of latent infectivity profile `profile`.
        void build_latent(const str_table &input, const InfecProfile &profile)
        {
            latent.clear();
            latent.append(input);
            latent.append(profile.list_latent_definitions(infec_profile));
        }
    };

    /// Memory chains run in, one at a time: model, modules built from definitions, totals, adaptation.
    /// * modules point into each other -> never copied, rebuilt.
    struct Workspace
    {
        /// Own `model_ptr`, register its keys in every key holder, bind `defs`.
        /// * model -> bound to this workspace, then its requirements declared;
        /// * `shared` -> its imported tables shared, read-cells unset until evaluation setup.
        Workspace(std::unique_ptr<ModelBase> model_ptr,
                  const Definitions &defs,
                  const Settings &settings,
                  const Data *shared = nullptr);
        Workspace(const Workspace &) = delete;
        Workspace &operator=(const Workspace &) = delete;
        ~Workspace();

        //////  Build  //////

        /// Build the modules that read slot `slot` of definitions.
        /// * import -> the slot is set first, then built on the first workspace;
        /// * rebuilt workspace -> every slot built again, in `Slot` order;
        /// * equivalence classes -> left to the caller, partitioned again at evaluation setup;
        /// * clears the evaluated flag.
        /// @note Every other workspace member a slot sets is set here, so imports and rebuilds agree.
        void build(Definitions::Slot slot);

        //////  Modules  //////

        /// Transmission model, bound to this workspace.
        std::unique_ptr<ModelBase> model;
        /// Definitions the modules are built from.
        const Definitions &defs;
        /// Imported tables, shared across engines, and run state.
        Data data;
        Parameter param;
        Rates rates;
        InfecProfile infec_profile;
        Augmented augm;
        Translator translator;

        //////  Totals and adaptation  //////

        /// Log-likelihood totals over households.
        LogLikTotal tot_log_lik;
        /// Route of each sampled parameter, by key.
        std::vector<ParamRoute> routes = {};
        /// Blocks of the run, see Rates::list_param_rw_blocks().
        std::vector<ParamBlock> pblocks = {};
        /// Whether each parameter belongs to a block.
        std::vector<char> in_pblock = {};
        /// Whether the adaptation warm-up runs, block members moving one at a time.
        bool is_warm = false;
        /// Robbins-Monro gain of the scalar proposals at this step.
        double gain = 0.;
        /// Robbins-Monro gain of the block factors at this step.
        double block_gain = 0.;

        /// Chain the workspace runs, index of its streams and of its posterior column block.
        Settings::Chain id = 0;

        //////  Setter  //////

        /// Initialize equivalence-class compression.
        /// * latent definitions required, not imported yet -> one class per household.
        void init_equivalence_class_compression();

        /// Throw if the generic profile transmits at generation time 0.
        /// * status-specific profile -> checked per member, with its variable.
        void assert_valid_infectivity_profile() const;

        //////  Getter  //////

        /// Latent `a` of member `m`.
        double aug(size_t a, const Individual &m) const { return augm.at(a, m).value(); }

        /// Data column names shared across simulations.
        std::vector<std::string> get_static_data_names() const;
        /// Column names of the per-replicate shallow list.
        std::vector<std::string> get_col_names() const;
        /// Per-household row offsets, so parallel simulations write disjoint rows.
        std::vector<size_t> get_row_offsets() const;

        //////  Log-likelihood  //////

        /// Observation-model log-likelihood of member `m`.
        const LogProba &obs_log_lik(const Individual &m) const noexcept { return augm.state.of(m).obs_log_lik; }
        /// Latent-model log-likelihood of member `m`.
        const LogProba &latent_log_lik(const Individual &m) const noexcept { return augm.state.of(m).latent_log_lik; }
        /// Sum the observation log-likelihood over infected members of `house`.
        double house_obs_log_lik(const Household &house) const
        {
            double obs = 0.;
            for (const Individual &indiv : house.infected_indivs())
                obs += obs_log_lik(indiv).value();
            return obs;
        }
        /// Sum the latent log-likelihood over infected members of `house`.
        double house_latent_log_lik(const Household &house) const
        {
            double model = 0.;
            for (const Individual &indiv : house.infected_indivs())
                model += latent_log_lik(indiv).value();
            return model;
        }
        /// Sum the transmission, observation and latent log-likelihood stored for `house`.
        /// * observation, latent -> only when the model declares them.
        double house_log_lik(const Household &house) const
        {
            double value = data.state.transm_log_lik(house).value();
            if (augm.scheme().is_any_obs_model())
                value += house_obs_log_lik(house);
            if (augm.scheme().is_any_latent_model())
                value += house_latent_log_lik(house);
            return value;
        }
    };
}

//////////////////////////////  MCMCBase  //////////////////////////////

namespace housepi
{
    class MCMCBase
    {
    public:
        Settings settings;
        Posterior posterior;
        /// Definitions of the run, read by every chain.
        Definitions defs;

        //////  Constructor  //////

        /// Build the first workspace on `model`.
        explicit MCMCBase(std::unique_ptr<ModelBase> model);
        virtual ~MCMCBase();

        //////  Getter  //////

        /// Workspace built by imports, read by evaluation and simulation.
        Workspace &chain() noexcept { return *workspaces_.front(); }
        const Workspace &chain() const noexcept { return *workspaces_.front(); }

    protected:
        /// Workspaces the chains run on, one chain at a time per workspace.
        /// * first -> built by imports, the only workspace outside a run;
        /// * others -> rebuilt from definitions for one run, see MCMC::prepare_inference().
        std::vector<std::unique_ptr<Workspace>> workspaces_ = {};
    };
}
