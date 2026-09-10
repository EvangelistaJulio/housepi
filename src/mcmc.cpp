/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/keys.h"
#include "tools/logproba.h"
#include "tools/messages.h"
#include "tools/nan.h"
#include "tools/parallel.h"
#include "tools/parser.h"
#include "tools/string.h"
#include "tools/types.h"
#include "tools/vector.h"

#include "augmented.h"
#include "contact_risk.h"
#include "data_wide.h"
#include "data_long.h"
#include "distribution.h"
#include "function.h"
#include "group_risk.h"
#include "household.h"
#include "individual.h"
#include "infec_profile.h"
#include "mcmc.h"
#include "mcmc_base.h"
#include "model.h"
#include "parameter.h"
#include "registry.h"
#include "setup.h"
#include "translator.h"

namespace housepi
{
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////  *** registry::RegModels : class ***  ////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace registry
    {
        RegModels &RegModels::instance()
        {
            static RegModels _inst;
            return _inst;
        }

        std::unique_ptr<ModelBase> RegModels::create(std::string_view name) const
        {
            auto it = factory_map_.find(name);
            if (it == factory_map_.end())
                msg::stop(
                    "Model registration failed: unknown model " + msg::val(std::string(name)) + ".\n");
            return it->second();
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  *** Translator : struct ***  //////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    const Attribute &Translator::to_attribute(
        const parse::Entry &node,
        const Individual *indiv_ptr) const
    {
        storage_ = Attribute(); // reset auxiliary storage to default empty

        if (!node.is_terminal())
        {
            helper_store_read_function_(node, indiv_ptr);
        }
        else if (parse::is_numeric(node.type))
        {
            storage_.init_numeric(str::tod(node.value, "argument with numeric type"));
        }
        else if (parse::is_parameter(node.type))
        {
            return param_.at(param_.khold(node.value));
        }
        else if (!node.is_empty())
        {
            if (!indiv_ptr)
            {
                std::string sms = "Individual specific definitions as " + node.err_str() +
                                  " are not supported here.";
                // add tip if latent (which does not require '@')
                if (parse::is_latent(node.type) && node.is_terminal())
                    sms += "\nDid you forget rlang injection operator " + msg::code("!!" + node.value) + "?";
                msg::stop(sms);
            }

            if (parse::is_data_wide(node.type))
            {
                // wire the column's read-cell: whatever produces the column points
                // that cell at itself, so the reader needs no rule of its own
                storage_.init_pointer(
                    &pdata_wide_->readcell(pdata_wide_->khold(node.value), *indiv_ptr));
                storage_.name(node.value, ErrType::Data);
            }
            else if (parse::is_latent(node.type))
            {
                return augm_.at(augm_.khold(node.value), *indiv_ptr);
            }
            else if (parse::is_group_risk(node.type))
            {
                return pgroup_risk_->at(pgroup_risk_->khold(node.value), *indiv_ptr);
            }
        }

        return storage_;
    }

    void Translator::helper_store_read_function_(
        const parse::Entry &node,
        const Individual *indiv_ptr) const
    {
        Attribute parent;
        parent.init_function(registry::create_function(
            node.func,
            ErrType::Function,
            {registry::Kind::Custom,
             registry::Kind::Transcendental,
             registry::Kind::Arithmetic,
             registry::Kind::Continuous,
             registry::Kind::Discrete,
             registry::Kind::InfecProfile}));

        for (const auto &child : node.children)
            parent.init_func_attr(child.key, to_attribute(child, indiv_ptr));

        storage_ = parent;
    }
    const Attribute &Translator::to_attribute(
        const std::string &expr,
        const Individual *indiv_ptr) const
    {
        return to_attribute(parse::Definition(expr).root, indiv_ptr);
    }

    const Attribute &Translator::to_numeric(
        const parse::Entry &node,
        const Individual *indiv_ptr) const
    {
        storage_ = to_attribute(node, indiv_ptr);
        if (!storage_.is_empty() && !storage_.is_constant())
            msg::stop(
                "Definition " + node.err_str() +
                " could not be interpreted as a numeric constant.");
        return storage_;
    }
    const Attribute &Translator::to_numeric(
        const std::string &expr,
        const Individual *indiv_ptr) const
    {
        return to_numeric(parse::Definition(expr).root, indiv_ptr);
    }

    void Translator::mark_param_produces_transm_cache(const parse::Entry &node) const
    {
        if (!parse::is_parameter(node.type))
            return;
        size_t p = param_.khold(node.value);
        if (!param_.at(p).is_constant())
            param_.mark_produces_transm_cache(p);
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  *** Posterior : struct ***  ///////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        // Help initializing parameter draws,
        void helper_init_param_draws(
            Posterior *post,
            const Parameter &param,
            unsigned long sample_size)
        {
            size_t n_sample = param.list_sample.size();
            size_t n_func = param.list_func.size();
            post->param_values.reserve(n_sample + n_func);
            post->param_priors.reserve(n_sample);
            post->param_names.reserve(n_sample + n_func);
            post->param_prior_names.reserve(n_sample);
            for (size_t p : param.list_sample)
            {
                post->param_values.emplace_back(sample_size, na::DOUBLE);
                post->param_priors.emplace_back(sample_size, na::DOUBLE);
                post->param_names.push_back(param.khold(p));
                post->param_prior_names.push_back("param_lprior_" + param.khold(p));
            }
            // Function-derived parameters: value trace only (no prior).
            for (size_t p : param.list_func)
            {
                post->param_values.emplace_back(sample_size, na::DOUBLE);
                post->param_names.push_back(param.khold(p));
            }
        }

        /// Help initializing household transmission log-likelihood draws.
        void helper_init_house_draws(
            Posterior *post,
            const DataWide &data_wide,
            unsigned long sample_size)
        {
            post->house_log_lik.reserve(data_wide.eq_class_size());
            for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
                post->house_log_lik.emplace_back(sample_size, na::DOUBLE);

            size_t n_houses = data_wide.report_house.total;
            post->house_names.reserve(n_houses);
            for (const Household &house : data_wide)
                post->house_names.push_back(data_wide.house_id(house));
        }

        // Help initializing latent draws,
        void helper_init_latent_draws(
            Posterior *post,
            const Augmented &augm,
            const DataWide &data_wide,
            const Setup &setup,
            unsigned long sample_size)
        {
            post->tracked_indivs = augm.init_individual_sample(
                setup.id_draws, setup.house_draws, setup.prng());

            if (post->tracked_indivs.empty())
                return;

            post->latent_values.resize(post->tracked_indivs.size());
            for (size_t i = 0; i < post->tracked_indivs.size(); ++i)
            {
                const Individual &indiv = *post->tracked_indivs[i];
                const auto &draws_list = augm.list_draws(indiv);
                auto &indiv_matrix = post->latent_values[i];
                indiv_matrix.reserve(draws_list.size());
                for (size_t a : draws_list)
                {
                    indiv_matrix.emplace_back(sample_size, na::DOUBLE);
                    post->latent_col_names.push_back(data_wide.id(indiv) + "_" + augm.khold(a));
                }
            }
        }
    }

    void Posterior::init(MCMCBase &mcmc)
    {
        count_ = 0;
        is_power_post_ = false;
        tot_log_lik.clear();
        param_values.clear();
        param_priors.clear();
        house_log_lik.clear();
        latent_values.clear();
        tracked_indivs.clear();
        param_names.clear();
        param_prior_names.clear();
        latent_col_names.clear();
        house_names.clear();

        const Setup &setup = mcmc.setup;
        unsigned long sample_size = setup.c_sample_size * setup.n_chains;
        is_power_post_ = setup.is_power_post();

        // power posteriors only track the joint log-likelihood total for the
        // thermodynamic ladder; ordinary inference never reads it, so it is the
        // single quantity stored in that mode (no pointwise / parameter draws).
        if (is_power_post_)
        {
            tot_log_lik.reserve(1);
            tot_log_lik.emplace_back(sample_size, na::DOUBLE);
            return;
        }

        helper_init_param_draws(this, mcmc.param, sample_size);
        helper_init_house_draws(this, mcmc.data_wide(), sample_size);
        helper_init_latent_draws(this, mcmc.augm, mcmc.data_wide(), setup, sample_size);
    }

    void Posterior::store(MCMCBase &mcmc)
    {
        // power posteriors store only the joint (tempered) log-likelihood that
        // feeds the thermodynamic ladder: transmission plus observation
        if (is_power_post_)
        {
            double value = mcmc.tot_transm_log_lik.value();
            if (mcmc.augm.is_any_obs_model())
                value += mcmc.tot_obs_log_lik.value();
            tot_log_lik[0][count_++] = value;
            return;
        }

        // parameter values and priors
        const Parameter &param = mcmc.param;
        size_t i = 0;
        for (size_t p : param.list_sample)
        {
            const Variable &var = param.at(p);
            param_values[i][count_] = var.value();
            param_priors[i++][count_] = var.log_pdf();
        }
        // function-derived parameters: store value, no prior
        for (size_t p : param.list_func)
            param_values[i++][count_] = param.at(p).value();

        // pointwise per-household log-likelihood:
        // transmission + observation + latent
        // NOTE: eq-classes internally set to distinct elsewhere, if const models
        const DataWide &data_wide = mcmc.data_wide();
        const bool with_obs = mcmc.augm.is_any_obs_model();
        const bool with_model = mcmc.augm.is_any_latent_model();
        for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
        {
            const Household &head = data_wide.eq_class_head(R);
            double value = data_wide.house_transm_log_lik(R);
            if (with_obs)
                value += mcmc.house_obs_log_lik(head);
            if (with_model)
                value += mcmc.house_latent_log_lik(head);
            house_log_lik[R][count_] = value;
        }

        // latent variable draws
        const Augmented &augm = mcmc.augm;
        for (size_t j = 0; j < tracked_indivs.size(); ++j)
        {
            const Individual &indiv = *tracked_indivs[j];
            auto &indiv_matrix = latent_values[j];
            size_t c = 0;
            for (size_t a : augm.list_draws(indiv))
                indiv_matrix[c++][count_] = mcmc.aug(a, indiv);
        }

        // update sample count
        ++count_;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////  *** MCMCBase : class ***  ////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    MCMCBase::MCMCBase()
        : pDataStore(std::make_shared<DataStore>()),
          augm(this), translator(param, augm) {}
    MCMCBase::~MCMCBase() = default;

    void MCMCBase::set_data_store(std::shared_ptr<DataStore> ds)
    {
        pDataStore = std::move(ds);
        translator.bind_data_wide(&pDataStore->data_wide);
    }

    void MCMCBase::init_key_holder(
        const std::vector<std::string> &data_wide_keys,
        const std::vector<std::string> &hh_member_keys,
        const std::vector<std::string> &hh_member_vec_keys,
        const std::vector<std::string> &hh_double_vec_keys,
        const std::vector<std::string> &hh_double_mtx_keys,
        const std::vector<std::string> &data_long_keys,
        const std::vector<std::string> &param_keys,
        const std::vector<std::string> &group_risk_keys,
        const std::vector<std::string> &latent_keys)
    {
        data_wide().init_key_holder(data_wide_keys);
        data_wide().init_cache_key_holders(hh_member_keys,
                                           hh_member_vec_keys,
                                           hh_double_vec_keys,
                                           hh_double_mtx_keys);
        data_long().init_key_holder(data_long_keys);
        param.init_key_holder(param_keys);
        group_risk.init_key_holder(group_risk_keys);
        augm.init_key_holder(latent_keys);

        translator.bind_data_wide(&data_wide());
        translator.bind_group_risk(&group_risk);
    }

    /// Wheter hosuehold members i and j are equivalent in terms of risks.
    /// Two members are in the same if group-risk attributes and
    /// every contact-risk involving one matches the corresponding one involving the
    /// other, both outgoing and incoming, across all other members.
    bool MCMCBase::is_same_risk_partition(const Household &house,
                                          MemberIdx i,
                                          MemberIdx j) const
    {
        if (i == j)
            return true;

        const Individual &indiv_i = house.at(i);
        const Individual &indiv_j = house.at(j);

        // check all group risks
        if (!group_risk.same_group(indiv_i, indiv_j))
            return false;

        // check all contact risks
        if (!contact_risk.same_contacts(house, i, j))
            return false;

        return true;
    }

    void MCMCBase::init_equivalence_class_compression()
    {
        // unique equivalence class if sampling latent,
        // latent models or per-member observation nodes
        bool force_unique = augm.is_unique_eq_class();
        // however, false could also mean definitions not provided:
        // check if required and in such case set unique preventively
        if (!force_unique && augm.list_latent_export.size() == 0)
            force_unique = augm.khold().size() > 0;

        data_wide().init_eq_classes(force_unique);
        data_wide().allocate_transm_log_lik();
    }

    namespace infprof
    {
        /// Throw if the profile transmits at generation time 0.
        /// UPDATE: unresolved argument per member leavs pdf undetermined
        /// only a determinate value is controlled
        void assert_zero_at_origin(const Function &f, bool is_discr_requested)
        {
            // skip for hand-written infectivity profiles
            if (!is_discr_requested && registry::instance(registry::Kind::InfecProfile).contains(f.name()))
                return;

            const Domain &domain = f.domain();
            // NOTE: not support(), whose dense rule would relax a floor at zero
            const Codomain reach = domain.relocated();
            // the bounds keep every attainable draw past the origin
            if (reach.is_positive())
                return;

            std::string sms = "Infectivity profile " + msg::val(f.name());

            if (na::is(domain.delocate(0.)))
            {
                // the origin reads an undrawn value the bounds do not keep away from it
                sms += " is not kept away from generation " +
                       std::string(is_discr_requested ? "day" : "time") + " 0.\n" +
                       "Its support is " + reach.str() + ", so it can put infectivity on the " +
                       (is_discr_requested ? "day" : "instant") + " of infection.";
            }
            else
            {
                const double y = f.pdf(0.);
                if (std::abs(y) < 1.e-10 || na::is(y))
                    return;

                sms += " is non-zero at generation ";
                if (is_discr_requested)
                    sms += "day 0 (F[0] = " + str::format(y) + ").\n" +
                           "A discrete-time profile admits no transmission on the day of infection.";
                else
                    sms += "time 0 (f(0) = " + str::format(y) + ").\n" +
                           "A continuous-time profile admits no transmission at the instant of infection.";
            }

            msg::stop(
                sms +
                "\nConsider truncating the distribution with " + msg::arg("min") +
                ", or shifting with " + msg::arg("loc") + ".");
        }
    }

    void MCMCBase::assert_valid_infectivity_profile() const
    {
        if (!infec_profile.generic && !infec_profile.is_latent())
            return;

        bool is_discr = infec_profile.is_required_discrete();

        if (infec_profile.generic)
            infprof::assert_zero_at_origin(*infec_profile.generic, is_discr);
        else
            for (const Household &house : data_wide())
                for (const Individual &indiv : house.infected_indivs())
                {
                    // check status-specific profiles
                    const Variable &avar = augm.at(infec_profile.latent_key, indiv);
                    // UPDATE: am checking it's a function somewhere?
                    infprof::assert_zero_at_origin(*avar.function(), is_discr);
                }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////  MCMCBase : monte carlo  /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    void MCMCBase::init_monte_carlo_()
    {
        // check if there is actually something to sample
        if (!param.is_any_sample() && !augm.is_any_sample())
            msg::stop(
                "Nothing to sample.\nCheck your parameter and latent definitions.");

        // whether adaptation is required
        bool is_adaptive = !param.list_adaptive.empty() || augm.is_any_adaptive();
        // initialize setup for mcmc
        setup.init_mcmc_configuration(is_adaptive);

        // set moving average smoothing factor for every Variable
        Variable::set_adapt_smoothing(1. / setup.n_adapt_window);
    }

    void MCMCBase::init_monte_carlo_chain_()
    {
        setup.prog_set_text();       // set progress bar text
        param.restart(setup.prng()); // reset parameters
        augm.restart(setup.prng());  // reset latent

        // re-initialize the adaptable proposal of every variable
        for (size_t p : param.list_adaptive)
            param.at(p).init_proposal();
        if (augm.is_any_adaptive())
            for (const Household &house : data_wide())
                for (const Individual &indiv : house.infected_indivs())
                    for (size_t a : augm.list_adaptive(indiv))
                        augm.at(a, indiv).init_proposal();
    }

    void MCMCBase::recover_total_log_lik_(bool with_cache,
                                          bool with_latent_ll,
                                          bool with_obs_ll)
    {
        tot_transm_log_lik.recover();
        if (with_obs_ll)
            tot_obs_log_lik.recover();
        if (with_latent_ll)
            tot_latent_log_lik.recover();

        // recover old log-likelihood values for households
        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
        {
            const Household &house = data_wide().eq_class_head(R);
            recover_transm_log_lik_(house);

            // recover from commit_total_log_posterior_
            if (with_latent_ll)
                for (const Individual &indiv : house.infected_indivs())
                    recover_latent_log_lik_(indiv);
            if (with_obs_ll)
                for (const Individual &indiv : house.infected_indivs())
                    recover_obs_log_lik_(indiv);
        }

        if (with_cache)
        {
            data_wide().recover_old_cached();
            infec_profile.recover();
        }
    }

    std::vector<std::string> MCMCBase::get_static_data_names() const
    {
        // exclude exported latent and observation-node columns: the latter are
        // generated by the simulation, so the generated value wins over static data.
        // surrogate columns are likewise produced (exported as columns of their own),
        // so they must not also surface as a static data column.
        std::vector<std::string> produced = vect::concatenate(
            augm.khold(augm.list_latent_export), augm.khold(augm.list_obs_member));
        for (const Augmented::Surrogate &surr : augm.list_surrogate)
            if (Augmented::is_own_column(surr))
                produced.push_back(augm.surrogate_name(surr));
        return vect::difference(data_wide().khold().all(), produced);
    }

    std::vector<size_t> MCMCBase::get_row_offsets() const
    {
        size_t n_houses = data_wide().report_house.total;
        std::vector<size_t> output(n_houses + 1, 0);
        for (size_t H = 0; H < n_houses; ++H)
            output[H + 1] = output[H] + data_wide().at(H).size();
        return output;
    }

    double MCMCBase::house_obs_log_lik(const Household &house) const
    {
        double obs = 0.;
        for (const Individual &indiv : house.infected_indivs())
            obs += obs_log_lik_(indiv).value();
        return obs;
    }

    double MCMCBase::house_latent_log_lik(const Household &house) const
    {
        double model = 0.;
        for (const Individual &indiv : house.infected_indivs())
            model += latent_log_lik_(indiv).value();
        return model;
    }

    std::vector<std::string> MCMCBase::get_col_names() const
    {
        std::vector<std::string> output = {"id", "house_id"};
        vect::append(output, get_static_data_names());
        vect::append(output, {"status"});

        output.reserve(output.size() +
                       augm.list_latent_export.size() +
                       augm.list_obs_member.size() +
                       augm.list_surrogate.size());

        for (auto idx : augm.list_latent_export)
            output.push_back(augm.khold(idx));

        for (auto idx : augm.list_obs_member)
            output.push_back(augm.khold(idx));

        for (const Augmented::Surrogate &surr : augm.list_surrogate)
            if (Augmented::is_own_column(surr))
                output.push_back(augm.surrogate_name(surr));

        return output;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////  *** MCMC : class ***  /////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////

    MCMC::MCMC(std::string_view model_name)
        : model(registry::RegModels::instance().create(model_name))
    {
        init_key_holder(model->key_Data(),
                        model->key_HHmember(),
                        model->key_HHmembervec(),
                        model->key_HHdoublevec(),
                        model->key_HHdoublemtx(),
                        model->key_Long(),
                        model->key_Param(),
                        model->key_GroupRisk(),
                        model->key_Latent());
        model->bind_solver(this);
        model->init_model_requirements();
    }

    void MCMC::detach_data_store()
    {
        if (pDataStore.use_count() <= 1)
            return;
        auto new_ds = std::make_shared<DataStore>();
        new_ds->data_wide.init_key_holder(model->key_Data());
        new_ds->data_wide.init_cache_key_holders(
            model->key_HHmember(),
            model->key_HHmembervec(),
            model->key_HHdoublevec(),
            model->key_HHdoublemtx());
        new_ds->data_long.init_key_holder(model->key_Long());
        set_data_store(std::move(new_ds));
    }

    //////////////////////////
    //////  EVALUATION  //////
    //////////////////////////

    double MCMC::evaluate_log_lik(std::string_view house_id)
    {
        setup.n_chains = 0; // using '0 chains' as a state variable
        init_for_evaluation_();
        param.draw_sample(setup.prng());
        augm.draw_sample(setup.prng());

        // evaluation at current state
        LogLikAccumulator acc;
        if (na::is(house_id))
        {
            // evaluation across all households
            commit_total_log_posterior_(true);
            acc.obs_ll = tot_obs_log_lik.value();
            acc.transm_ll = tot_transm_log_lik.value();
            acc.latent_ll = tot_latent_log_lik.value();
        }
        else
        {
            // single household evaluation
            infec_profile.update_cache();
            acc = commit_house_log_posterior_(
                data_wide().at(house_id),
                true, augm.is_any_latent_model(), augm.is_any_obs_model());
        }

        // complete-data household log-likelihood (joint over latent):
        // observation + transmission + latent model
        return acc.obs_ll + acc.transm_ll + acc.latent_ll;
    }

    void MCMC::init_for_evaluation_()
    {
        DataWide &data_wide = this->data_wide();
        data_wide.assert_not_all_missing_status();

        // materialize surrogate data columns before anything reads the data store
        data_wide.reset_readcells(augm.get_produced_cols());
        augm.init_for_evaluation();
        data_wide.init_eq_classes(augm.is_unique_eq_class());
        data_wide.allocate_all_cache();
        data_wide.allocate_transm_log_lik();

        param.finalize();
        group_risk.finalize();
        contact_risk.finalize();
        augm.finalize();
        infec_profile.finalize();
        data_wide.finalize();

        param.assert_risk_identifiability(group_risk, contact_risk);
        assert_valid_infectivity_profile();

        model->init_model_for_evaluation();
        for (size_t R = 0; R < data_wide.eq_class_size(); ++R)
            model->init_house_for_evaluation(data_wide.eq_class_head(R));
    }

    void MCMC::commit_total_log_posterior_(bool with_cache,
                                           bool with_latent_ll,
                                           bool with_obs_ll)
    {
        DataWide &data_wide = this->data_wide();

        if (with_cache)
        {
            data_wide.save_all_cached();
            infec_profile.update_cache();
        }

        // NOTE: eq-class coincide with all households when a model or an observation exist
        const bool do_latent = with_latent_ll && augm.is_any_latent_model();
        const bool do_obs = with_obs_ll && augm.is_any_obs_model();
        LogLikAccumulator acc = parallel::for_reduce<LogLikAccumulator>(
            data_wide.eq_class_size(),
            [&](size_t R) -> LogLikAccumulator
            {
                double w = data_wide.eq_class_weight(R);
                LogLikAccumulator a = commit_house_log_posterior_(
                    data_wide.eq_class_head(R),
                    with_cache, do_latent, do_obs);
                return {w * a.obs_ll, w * a.transm_ll, w * a.latent_ll};
            });

        tot_transm_log_lik.set(acc.transm_ll);
        if (do_obs)
            tot_obs_log_lik.set(acc.obs_ll);
        if (do_latent)
            tot_latent_log_lik.set(acc.latent_ll);
    }

    LogLikAccumulator MCMC::commit_house_log_posterior_(const Household &house,
                                                        bool with_cache,
                                                        bool with_latent_ll,
                                                        bool with_obs_ll)
    {
        LogLikAccumulator acc;

        // re-compute transmission cached quantities before the commit
        if (with_cache)
            model->update_house_transmission_cache(house);

        // --- Transmission ---
        acc.transm_ll = commit_transm_log_lik_(house);
        // --- Observation ---
        if (with_obs_ll)
            for (const Individual &indiv : house.infected_indivs())
                acc.obs_ll += augm.commit_obs_log_lik(indiv);
        // --- Latent ---
        if (with_latent_ll)
            for (const Individual &indiv : house.infected_indivs())
                acc.latent_ll += augm.commit_latent_log_lik(indiv);
        return acc;
    }

    double MCMC::commit_total_latent_log_lik_ratio_()
    {
        double acc = parallel::for_reduce<double>(
            data_wide().eq_class_size(),
            [&](size_t R) -> double
            {
                const Household &house = data_wide().eq_class_head(R);
                double amodel = 0.;
                for (const Individual &indiv : house.infected_indivs())
                    amodel += augm.commit_latent_log_lik(indiv);
                return data_wide().eq_class_weight(R) * amodel;
            });
        tot_latent_log_lik.set(acc);
        return tot_latent_log_lik.ratio();
    }

    void MCMC::recover_total_latent_log_lik_()
    {
        tot_latent_log_lik.recover();
        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
        {
            const Household &house = data_wide().eq_class_head(R);
            for (const Individual &indiv : house.infected_indivs())
                recover_latent_log_lik_(indiv);
        }
    }

    /////////////////////////////
    ///  METROPOLIS-HASTINGS  ///
    /////////////////////////////

    void MCMC::run_inference()
    {
        init_monte_carlo_();
        init_for_evaluation_();
        setup.prog_open_for_mcmc_sample();
        run_monte_carlo_();
        setup.prog_close();
    }

    void MCMC::run_monte_carlo_()
    {
        setup.prog_set_text(); // refresh label
        posterior.init(*this); // initialize posterior sample

        // run multiple MCMC chains
        for (setup.chain_id = 0; setup.chain_id < setup.n_chains; ++setup.chain_id)
        {
            init_monte_carlo_chain_();
            init_guarded_commit_total_();
            if (setup.inv_temp == 0.0)
            {
                run_mc_finite_adaptation_<true>();
                run_mc_posterior_sampling_<true>();
            }
            else
            {
                run_mc_finite_adaptation_<false>();
                run_mc_posterior_sampling_<false>();
            }
        }
    }

    void MCMC::init_guarded_commit_total_()
    {
        data_wide().save_all_cached();
        infec_profile.update_cache();

        // single serial pass
        LogLikAccumulator acc;
        for (size_t R = 0; R < data_wide().eq_class_size(); ++R)
        {
            double w = data_wide().eq_class_weight(R);
            LogLikAccumulator a = init_guarded_commit_house_(data_wide().eq_class_head(R));
            acc += {w * a.obs_ll, w * a.transm_ll, w * a.latent_ll};
        }

        tot_transm_log_lik.set(acc.transm_ll);
        if (augm.is_any_obs_model())
            tot_obs_log_lik.set(acc.obs_ll);
        if (augm.is_any_latent_model())
            tot_latent_log_lik.set(acc.latent_ll);
    }

    LogLikAccumulator MCMC::init_guarded_commit_house_(Household &house)
    {
        LogLikAccumulator acc;

        // --- Transmission ---
        model->update_house_transmission_cache(house);
        acc.transm_ll = commit_transm_log_lik_(house);
        if (data_wide().transm_log_lik(house).is_inf())
            msg::stop(
                "Initial transmission log-likelihood returned infinity for " +
                data_wide().str(house) + ".");

        // --- Observation ---
        // resample each individual variable upon failure
        if (augm.is_any_obs_model())
            for (const Individual &indiv : house.infected_indivs())
            {
                if (augm.list_obs_models(indiv).empty())
                    continue;

                size_t attempt = 0;
                while (!std::isfinite(augm.eval_log_obs(indiv)))
                {
                    ++attempt;
                    if (!augm.is_any_sample() || attempt >= N_INIT_ATTEMPTS_)
                        augm.throw_obs_init_error(indiv);
                    augm.restart(indiv, setup.prng());
                }
                acc.obs_ll += augm.commit_obs_log_lik(indiv);
            }

        // --- Latent ---
        if (augm.is_any_latent_model())
            for (const Individual &indiv : house.infected_indivs())
            {
                if (augm.list_latent_models(indiv).empty())
                    continue;

                acc.latent_ll += augm.commit_latent_log_lik(indiv);
                if (latent_log_lik_(indiv).is_inf())
                    msg::stop(
                        "Initial latent log-likelihood returned infinity for " +
                        data_wide().str(indiv) + ".");
            }

        return acc;
    }

    template <bool Beta0>
    void MCMC::run_mc_finite_adaptation_()
    {
        for (unsigned long i = 0; i < setup.n_adapt; ++i)
        {
            move_all_<true, Beta0>();

            // vanishing ratio of the adaptation budget already spent
            const double ratio = static_cast<double>(i) / setup.n_adapt;

            for (size_t p : param.list_adaptive)
                param.at(p).adapt_proposal(ratio);
            if (augm.is_any_adaptive())
                for (const Household &house : data_wide())
                    for (const Individual &indiv : house.infected_indivs())
                        for (size_t a : augm.list_adaptive(indiv))
                            augm.at(a, indiv).adapt_proposal(ratio);
        }
    }

    template <bool Beta0>
    void MCMC::run_mc_posterior_sampling_()
    {
        for (unsigned long i = 0; i < setup.n_burn; ++i)
            move_all_<false, Beta0>();

        for (unsigned long i = setup.n_burn; i < setup.n_steps; i += setup.n_thin)
        {
            for (unsigned long t = 1; t <= setup.n_thin; ++t)
                move_all_<false, Beta0>();
            // NOTE: recomputing likelihood, not tracked at inverse temperature 0
            if constexpr (Beta0)
                commit_total_log_posterior_(true, false, augm.is_any_obs_model());
            posterior.store(*this);
        }
    }

    template <bool Adapt, bool Beta0>
    void MCMC::move_all_()
    {
        // move latent variables to sample, if any
        move_all_latent_<Adapt, Beta0>();

        // move standard model parameters
        for (size_t p : param.list_sample_standard)
            move_single_parameter_<Adapt, Beta0>(p, false);
        // move parameters with cached dependencies
        for (size_t p : param.list_sample_cached)
            move_single_parameter_<Adapt, Beta0>(p, true);

        // update progress bar
        setup.prog_update();
    }

    template <bool Adapt, bool Beta0>
    void MCMC::move_single_parameter_(size_t p, bool with_cache)
    {
        Variable &this_par = param.at(p);
        // set new value of the parameter and store the old to buffer
        this_par.propose(setup.prng());

        // proceed only if the new value is consistent
        // otherwise automatically refuse move
        if (this_par.is_valid_move())
        {
            const bool with_latent_ll = augm.is_any_latent_model() &&
                                        augm.param_produces_latent_log_lik(p);

            /////////////////////////////
            ///  METROPOLIS-HASTINGS  ///
            /////////////////////////////

            if constexpr (Beta0)
            {
                // at inverse temperature 0 the likelihood totals are not tracked
                // they are refreshed from scratch before each stored draw
                double r = this_par.get_move_log_ratio();
                if (with_latent_ll)
                    r += commit_total_latent_log_lik_ratio_();
                if (setup.exponential() > -r)
                {
                    if constexpr (Adapt)
                        this_par.commit_proposal_move();
                    return;
                }
                if (with_latent_ll)
                    recover_total_latent_log_lik_();
            }
            else
            {
                const bool with_obs_ll = augm.is_any_obs_model() &&
                                         augm.param_produces_obs_log_lik(p);

                commit_total_log_posterior_(with_cache, with_latent_ll, with_obs_ll);

                // (tempered) Metropolis-Hastings acceptance log-ratio
                double r = tot_transm_log_lik.ratio();
                if (with_obs_ll)
                    r += tot_obs_log_lik.ratio();
                r *= setup.inv_temp;
                if (with_latent_ll)
                    r += tot_latent_log_lik.ratio();
                r += this_par.get_move_log_ratio();

                if (setup.exponential() > -r)
                {
                    if constexpr (Adapt)
                        this_par.commit_proposal_move();
                    return;
                }

                // else move is refused and log-likelihood recovered
                recover_total_log_lik_(with_cache, with_latent_ll, with_obs_ll);
            }
        }
        this_par.restore();
        if constexpr (Adapt)
            this_par.refuse_proposal_move();
    }

    template <bool Adapt, bool Beta0>
    void MCMC::move_all_latent_()
    {
        if (!augm.is_any_sample())
            return;

        // run per-house latent moves, and accumulate accepted log-density deltas
        LogLikAccumulator acc = parallel::for_reduce<LogLikAccumulator>(
            data_wide().report_house.total,
            [&](size_t H) -> LogLikAccumulator
            { return move_latent_single_house_<Adapt, Beta0>(H); });
        // at inverse temperature 0 the likelihood totals are not tracked
        // they are refreshed from scratch before each stored draw
        if constexpr (!Beta0)
        {
            if (augm.is_any_obs_model())
                tot_obs_log_lik.increment(acc.obs_ll);
            tot_transm_log_lik.increment(acc.transm_ll);
        }
        tot_latent_log_lik.increment(acc.latent_ll);
    }

    template <bool Adapt, bool Beta0>
    LogLikAccumulator MCMC::move_latent_single_house_(size_t H)
    {
        const Household &house = data_wide().at(H);

        LogLikAccumulator acc;
        for (const Individual &indiv : house.infected_indivs())
            for (size_t a : augm.list_sample(indiv))
                acc += move_latent_single_variable_<Adapt, Beta0>(house, indiv, a);
        return acc;
    }

    template <bool Adapt, bool Beta0>
    LogLikAccumulator MCMC::move_latent_single_variable_(const Household &house,
                                                         const Individual &indiv,
                                                         size_t a)
    {
        Variable &avar = augm.at(a, indiv);

        // keep track of model-consumers
        // NOTE: if 'a' is producer of other latent models, their model-ratio may change
        const bool has_consumers = augm.has_model_consumers(indiv, a);
        double cons_ll_old = 0.;
        if (has_consumers)
            cons_ll_old = augm.sum_guarded_latent_log_lik(indiv, a);

        // whether latent variable is a producer of transmision (and/or its cache)
        const bool produces_transm = augm.produces_transm(indiv, a);

        avar.propose(setup.prng()); // propose new value and save old
        if (avar.is_valid_move())   // proceed only if new value is valid
        {
            // recompute the cached quantities consuming this member
            if (produces_transm)
            {
                data_wide().save_member_cache(house, indiv);
                model->update_member_transmission_cache(house, indiv);
            }

            /////////////////////////////
            ///  METROPOLIS-HASTINGS  ///
            /////////////////////////////

            const double d_latent = has_consumers
                                        ? augm.commit_latent_log_lik_ratio_with_consumers(indiv, a, cons_ll_old)
                                        : commit_latent_log_lik_ratio(indiv, a);

            // NOTE: likelihood terms are tempered, the latent model is not.
            if constexpr (Beta0)
            {
                // at inverse temperature 0 the likelihood totals are not tracked
                // they are refreshed from scratch before each stored draw
                if (setup.exponential() > avar.proposal_log_ratio() - d_latent)
                {
                    if constexpr (Adapt)
                        avar.commit_proposal_move();
                    return {0., 0., d_latent};
                }
            }
            else
            {
                const double d_obs = commit_obs_log_lik_ratio(indiv, a);
                const double d_transm = produces_transm ? commit_transm_log_lik_ratio_(house) : 0.;

                const double r = setup.inv_temp * (d_obs + d_transm) + d_latent - avar.proposal_log_ratio();
                if (setup.exponential() > -r)
                {
                    if constexpr (Adapt)
                        avar.commit_proposal_move();
                    return {d_obs, d_transm, d_latent};
                }

                if (augm.produces_obs(indiv, a))
                    recover_obs_log_lik_(indiv);
                if (produces_transm)
                    recover_transm_log_lik_(house);
            }
            recover_latent_log_lik_(indiv);

            if (produces_transm)
                data_wide().recover_member_old_cached(house, indiv);
        }
        avar.restore();
        if constexpr (Adapt)
            avar.refuse_proposal_move();

        return {0., 0., 0.};
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////  MCMC : estimate log_evidence  ////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    namespace
    {
        /// Welford online's algorithm for mean and variance.
        /// Based on @cite welford_note_1962
        /// A single non-finite draw makes the mean infinite and variance undefined.
        std::pair<double, double> compute_stats(const std::vector<double> &data)
        {
            size_t n = 0;
            double avg = 0., sqr_diff = 0.;
            for (double v : data)
            {
                if (!std::isfinite(v))
                    return {-inf::DOUBLE, na::DOUBLE};

                const double delta = v - avg;
                avg += delta / ++n;
                sqr_diff += delta * (v - avg);
            }

            return {avg, n > 1 ? sqr_diff / static_cast<double>(n - 1) : 0.};
        }
    }

    PowerLadder::PowerLadder(size_t n_rungs)
    {
        size_t n = n_rungs + 1;
        inv_temp.reserve(n);
        mean_log_lik.reserve(n);
        var_log_lik.reserve(n);
        draws.reserve(n);
    }

    void PowerLadder::push(double inv_t, std::vector<double> row)
    {
        const auto [avg, var] = compute_stats(row);
        inv_temp.push_back(inv_t);
        mean_log_lik.push_back(avg);
        var_log_lik.push_back(var);
        draws.push_back(std::move(row));
    }

    void PowerLadder::insert(size_t i, double inv_t, std::vector<double> row)
    {
        const auto [avg, var] = compute_stats(row);
        inv_temp.insert(inv_temp.begin() + i, inv_t);
        mean_log_lik.insert(mean_log_lik.begin() + i, avg);
        var_log_lik.insert(var_log_lik.begin() + i, var);
        draws.insert(draws.begin() + i, std::move(row));
    }

    PowerLadder MCMC::run_power_posteriors(std::string_view scheme, unsigned long n_rungs)
    {
        // check number of rungs for the thermodynamic ladder
        if (n_rungs < 2)
            msg::stop(
                "N. of rungs was not greater than 1.");

        setup.n_rungs = n_rungs;          // NOTE: set before progress bar
        init_monte_carlo_();              // initialize mcmc sampling
        init_for_evaluation_();           // initialize for evaluation
        setup.prog_open_for_power_post(); // print progress bar
        // UPDATE: always under estimate progress-bar time if starting from itm = 0

        //////////////////////////
        ///  POWER-POSTERIORS  ///
        //////////////////////////

        PowerLadder ladder(n_rungs);
        power_post_zero_(ladder);

        // NOTE: adaptive is undefined if -Int at inverse-temperature = 0, fall back
        ladder.scheme = scheme;
        if (scheme == "adaptive" && !std::isfinite(ladder.mean_log_lik.front()))
        {
            msg::warn(
                "Adaptive scheme is undefined at inverse temperature 0 with -Inf log-likelihood.\n"
                "Falling back to the power-fractions scheme.");
            ladder.scheme = "power-fractions";
        }

        if (ladder.scheme == "adaptive")
            power_post_adaptive_(ladder);
        else
            power_post_fractions_(ladder);

        setup.prog_close(); // close progress bar

        // restore setup state variables
        setup.inv_temp = 1.0;
        setup.rung_id = 0;
        setup.n_rungs = 0;

        return ladder;
    }

    void MCMC::power_post_zero_(PowerLadder &ladder)
    {
        setup.rung_id = 0;
        setup.inv_temp = 0.0;
        run_monte_carlo_();
        ladder.push(setup.inv_temp, posterior.tot_log_lik[0]);
    }

    void MCMC::power_post_adaptive_(PowerLadder &ladder)
    {
        // inverse temperature 1
        setup.rung_id = 1;
        setup.inv_temp = 1.0;
        run_monte_carlo_();
        ladder.push(setup.inv_temp, posterior.tot_log_lik[0]);

        // symmetrized sum_i DKL between powerpost[i] and [i+1]
        // equals to sum_i Delta_inv_temp_i * (mean_log_lik[i+1] - mean_log_lik[i])
        // subdivide interval with the 'worst' contribution
        while (setup.rung_id < setup.n_rungs)
        {
            ++setup.rung_id;

            const auto &itm = ladder.inv_temp;
            const auto &avg = ladder.mean_log_lik;
            const auto &var = ladder.var_log_lik;

            size_t k = 0;
            double worst = -1.;
            for (size_t i = 0; i + 1 < itm.size(); ++i)
            {
                const double s = std::fabs((itm[i + 1] - itm[i]) * (avg[i + 1] - avg[i]));
                if (s > worst)
                {
                    worst = s;
                    k = i;
                }
            }

            // place new inverse temperature within [itm[k], itm[k + 1]]
            double new_t;
            if (!std::isfinite(avg[k]) || !std::isfinite(avg[k + 1]) ||
                !std::isfinite(var[k]) || !std::isfinite(var[k + 1]) ||
                avg[k + 1] <= avg[k])
                // midpoint if non-monotone (Monte Carlo error)
                // or non-finite statistic (e.g. itm = 0)
                // tangent formula is undefined
                new_t = 0.5 * (itm[k] + itm[k + 1]);
            else
            {
                // intersection of tangents at interval ends
                new_t = (avg[k + 1] - itm[k + 1] * var[k + 1] - avg[k] + itm[k] * var[k]) / (var[k] - var[k + 1]);
                // variance-weighted average if intersection falls outside interval
                if (new_t <= itm[k] || new_t >= itm[k + 1])
                    new_t = itm[k] + var[k + 1] * (itm[k + 1] - itm[k]) / (var[k] + var[k + 1]);
            }

            // power posterior sampling
            setup.inv_temp = new_t;
            run_monte_carlo_();
            ladder.insert(k + 1, setup.inv_temp, posterior.tot_log_lik[0]);
        }
    }

    void MCMC::power_post_fractions_(PowerLadder &ladder)
    {
        for (setup.rung_id = 1; setup.rung_id <= setup.n_rungs; ++setup.rung_id)
        {
            setup.inv_temp = std::pow(static_cast<double>(setup.rung_id) / static_cast<double>(setup.n_rungs), 5.);
            run_monte_carlo_();
            ladder.push(setup.inv_temp, posterior.tot_log_lik[0]);
        }
    }
}
