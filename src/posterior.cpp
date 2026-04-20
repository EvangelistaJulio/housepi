/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/nan.h"

#include "augmented.h"
#include "data_wide.h"
#include "household.h"
#include "individual.h"
#include "mcmc_base.h"
#include "parameter.h"
#include "setup.h"

namespace housepi
{
    // Help initializing parameter draws,
    void helper_init_param_draws(
        Posterior *post,
        const Parameter &param,
        unsigned long sample_size)
    {
        size_t n = param.list_sample.size();
        post->param_values.reserve(n);
        post->param_priors.reserve(n);
        post->param_names.reserve(n);
        post->param_prior_names.reserve(n);
        for (size_t p : param.list_sample)
        {
            post->param_values.emplace_back(sample_size, na::DOUBLE);
            post->param_priors.emplace_back(sample_size, na::DOUBLE);
            post->param_names.push_back(param.khold().as_string(p));
            post->param_prior_names.push_back("lprior_" + param.khold().as_string(p));
        }
    }

    /// Help initializing household log-likelihood draws.
    void helper_init_house_draws(
        Posterior *post,
        const DataWide &data_wide,
        unsigned long sample_size)
    {
        size_t n_restricted = data_wide.distinct_size();
        post->house_log_lik.reserve(n_restricted);
        for (size_t R = 0; R < n_restricted; ++R)
            post->house_log_lik.emplace_back(sample_size, na::DOUBLE);

        size_t n_houses = data_wide.report_house.total;
        post->house_names.reserve(n_houses);
        for (Household house : data_wide)
            post->house_names.push_back(house.house_id());
    }

    // Help initializing augmented draws,
    void helper_init_augm_draws(
        Posterior *post,
        const Augmented &augm,
        const std::vector<std::string> &id_draws,
        unsigned long sample_size)
    {
        post->augm_individuals = augm.init_individual_sample(id_draws);

        if (post->augm_individuals.empty())
            return;

        post->augm_values.resize(post->augm_individuals.size());
        for (size_t i = 0; i < post->augm_individuals.size(); ++i)
        {
            const Individual &member = *post->augm_individuals[i];
            const auto &draws_list = augm.list_draws(member);
            auto &member_matrix = post->augm_values[i];
            member_matrix.reserve(draws_list.size());
            for (size_t a : draws_list)
            {
                member_matrix.emplace_back(sample_size, na::DOUBLE);
                post->augm_col_names.push_back(member.id() + "_" + augm.khold().as_string(a));
            }
        }
    }

    void Posterior::init(MCMCBase &mcmc)
    {
        clear();

        const Setup &setup = mcmc.setup;
        unsigned long sample_size = setup.c_sample_size * setup.n_chains;
        is_power_post_ = setup.is_power_post();

        // always allocate total log-likelihood
        total_log_lik.reserve(1);
        total_log_lik.emplace_back(sample_size, na::DOUBLE);

        //// UPDATE: I would like to store also power posterior sample.
        /// NOTE: it will depend on chosen rungs in that case.
        if (is_power_post_)
            return;

        helper_init_param_draws(this, mcmc.param, sample_size);
        helper_init_house_draws(this, mcmc.data_wide(), sample_size);
        helper_init_augm_draws(this, mcmc.augm, setup.id_draws, sample_size);
    }

    void Posterior::store(MCMCBase &mcmc)
    {
        if (!is_power_post_)
        {
            // parameter values and priors
            const Parameter &param = mcmc.param;
            size_t i = 0;
            for (size_t p : param.list_sample)
            {
                const Variable &v = param.at(p);
                param_values[i][count_] = v.value();
                param_priors[i][count_] = v.get_distr_log_pdf();
                ++i;
            }

            // household log-likelihoods
            const DataWide &data_wide = mcmc.data_wide();
            for (size_t R = 0; R < data_wide.distinct_size(); ++R)
                house_log_lik[R][count_] = data_wide.house_log_lik(R);

            // augmented variable draws
            const Augmented &augm = mcmc.augm;
            for (size_t j = 0; j < augm_individuals.size(); ++j)
            {
                const Individual &member = *augm_individuals[j];
                auto &member_matrix = augm_values[j];
                size_t c = 0;
                for (size_t a : augm.list_draws(member))
                    member_matrix[c++][count_] = member.aug(a);
            }
        }

        // total log-likelihood (always stored, including power posterior)
        total_log_lik[0][count_] = mcmc.tot_log_lik.value();
        ++count_;
    }

    void Posterior::clear()
    {
        count_ = 0;
        is_power_post_ = false;

        for (auto &v : total_log_lik)
            std::vector<double>().swap(v);
        std::vector<std::vector<double>>().swap(total_log_lik);

        for (auto &v : param_values)
            std::vector<double>().swap(v);
        std::vector<std::vector<double>>().swap(param_values);

        for (auto &v : param_priors)
            std::vector<double>().swap(v);
        std::vector<std::vector<double>>().swap(param_priors);

        for (auto &v : house_log_lik)
            std::vector<double>().swap(v);
        std::vector<std::vector<double>>().swap(house_log_lik);

        for (auto &m : augm_values)
        {
            for (auto &v : m)
                std::vector<double>().swap(v);
            std::vector<std::vector<double>>().swap(m);
        }
        std::vector<std::vector<std::vector<double>>>().swap(augm_values);

        std::vector<const Individual *>().swap(augm_individuals);
        std::vector<std::string>().swap(param_names);
        std::vector<std::string>().swap(param_prior_names);
        std::vector<std::string>().swap(augm_col_names);
        std::vector<std::string>().swap(house_names);
    }
}
