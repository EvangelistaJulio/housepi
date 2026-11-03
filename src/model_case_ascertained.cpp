/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

//////////////////////////////  case_ascertained  //////////////////////////////

#include "housepi.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE case_ascertained
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()

/// Model ascertained cases with augmented infection times and a single introduction per household.
/// Simulation infects a member drawn uniformly at the time origin, see ModelBase::set_introduction().
/// Based on @cite cauchemez_bayesian_2004
class Model final : public ModelBase
{
public:
    HOUSEPI_DECLARE_LATENT_KEYS(infec_time)
    HOUSEPI_DECLARE_DATA_KEYS(end_date)
    HOUSEPI_DECLARE_GROUP_RATE_KEYS(community, susceptibility, infectivity)
    HOUSEPI_DECLARE_HOUSE_MEMBER_KEYS(intro_case)
    HOUSEPI_DECLARE_HOUSE_DOUBLE_MTX_KEYS(inst, cumul)
    HOUSEPI_DECLARE_PARAM_KEYS(alpha, beta, gamma, hh_size_ref)

    /// Declare data, parameter, rate and latent requirements.
    void init_model_requirements() override final
    {
        request_data_discrete(
            {Data::end_date});
        request_data_every_member(
            {Data::end_date});
        request_param_positive(
            {Param::alpha, Param::hh_size_ref});
        request_param_rw_block(
            {Param::beta, Param::gamma});
        request_rate_identifiability(
            Param::alpha, {GroupRate::community});
        request_rate_identifiability(
            Param::beta, {GroupRate::susceptibility, GroupRate::infectivity}, true);
        request_rate_status_stratified(
            {GroupRate::infectivity});
        request_latent_transmission(
            {Latent::infec_time});
#ifdef _HOUSEPI_DISCRETE_TIME
        request_latent_discrete(
            {Latent::infec_time});
        request_infec_profile_discrete();
#endif
    }

    /// Whether the model uses contact rates.
    bool has_contact_rate() const override final { return true; }

    /// Whether the model requires an infectivity profile definition.
    bool has_infec_profile() const override final { return true; }

    /// Whether households enter the study through a case.
    bool has_case_ascertainment() const override final { return true; }

    /// Assert no member has missing status.
    void init_model_for_evaluation() override final { assert_no_missing_status(); }

    /// Assert `house` has an infected member.
    void init_house_for_evaluation(const Household &house) override final
    {
        if (house.n_infected_status() == 0)
            msg::stop(
                "No infected members were found in " + str(house) +
                ".\nCase-ascertained model assumes at least one infected member per household.");
    }

    //////////////////////////////  Model : household log-likelihood  //////////////////////////////

    /// Return transmission rate, scaled by household size.
    /// Based on @cite cauchemez_bayesian_2004
    double transmission_rate(MemberIdx house_size) const
    {
        return par(Param::beta) *
               std::exp(par(Param::gamma) * std::log(par(Param::hh_size_ref) / house_size));
    }

    /// Return end of follow-up of `indiv`: `end_date` is the last fully observed day, included.
    double follow_up_end(const Individual &indiv) const
    {
#ifdef _HOUSEPI_DISCRETE_TIME
        return data(Data::end_date, indiv);
#else
        return data(Data::end_date, indiv) + 1.;
#endif
    }

    /// Compute transmission log-likelihood within household.
    double compute_house_transm_log_lik(const Household &house) override final
    {
        const MemberIdx i0 = cache_member(house, HouseMember::intro_case);
        const Individual &indiv_0 = house.at(i0);
        const double t_intro = aug(Latent::infec_time, indiv_0);

        const auto &inst_mtx = cache_double_matrix(house, HouseDoubleMtx::inst);
        const auto &cumul_mtx = cache_double_matrix(house, HouseDoubleMtx::cumul);

        const double alpha = par(Param::alpha);
        const double beta = transmission_rate(house.size());

        double output = 0.;
        for (const Individual &indiv_s : house.escaped_indivs())
        {
            const double exposure = follow_up_end(indiv_s) - t_intro;
            output -= alpha * exposure * group_rate(GroupRate::community, indiv_s);
        }

        for (MemberIdx i : house.infected_members())
        {
            const Individual &indiv_i = house.at(i);
            const double mu_i = group_rate(GroupRate::infectivity, indiv_i);

            for (MemberIdx s : house.escaped_members())
                output -= mu_i * beta * cumul_mtx[i][s] *
                          contact_rate(house, i, s) *
                          group_rate(GroupRate::susceptibility, house.at(s));

            if (i == i0)
                continue;

            // rate of infection from the other infected members j
            double roi_inst_i = 0.;
            double roi_cumul_i = 0.;
            for (MemberIdx j : house.infected_members())
            {
                if (j == i)
                    continue;

                const double mu_ji = group_rate(GroupRate::infectivity, house.at(j)) *
                                     contact_rate(house, j, i);
                roi_inst_i += mu_ji * inst_mtx[j][i];
                roi_cumul_i += mu_ji * cumul_mtx[j][i];
            }

            const double a_i = alpha * group_rate(GroupRate::community, indiv_i);
            const double b_i = beta * group_rate(GroupRate::susceptibility, indiv_i);

            roi_inst_i *= b_i;
            roi_cumul_i *= b_i;
            roi_inst_i += a_i;

#ifdef _HOUSEPI_DISCRETE_TIME
            // days escaped before infection: none for a member infected the day after the introduction
            // a co-primary cannot survive fewer
            const double survived = aug(Latent::infec_time, indiv_i) - 1. - t_intro;
            if (survived > 0.)
                roi_cumul_i += a_i * survived;
            output += math::log1mexp(roi_inst_i) - roi_cumul_i;
#else
            roi_cumul_i += a_i * (aug(Latent::infec_time, indiv_i) - t_intro);
            output += std::log(roi_inst_i) - roi_cumul_i;
#endif
        }

        return output;
    }

    /// Export community and person-to-person infection rates of `house`, by member.
    /// * `community_rate`: one per member;
    /// * `p2p_rate`: infector by row, infectee by column, diagonal 0, scaling the infector's profile.
    void export_house_quantities(const Household &house, HouseQuantities &out) override final
    {
        const MemberIdx house_size = house.size();
        const double alpha = par(Param::alpha);
        const double beta = transmission_rate(house_size);

        std::vector<double> commu(house_size);
        matrix p2p(house_size, std::vector<double>(house_size, 0.));
        for (MemberIdx i = 0; i < house_size; ++i)
        {
            const Individual &indiv_i = house.at(i);
            commu[i] = alpha * group_rate(GroupRate::community, indiv_i);
            const double inf_i = beta * group_rate(GroupRate::infectivity, indiv_i);
            for (MemberIdx j = 0; j < house_size; ++j)
                if (j != i)
                    p2p[i][j] = inf_i * contact_rate(house, i, j) * group_rate(GroupRate::susceptibility, house.at(j));
        }
        out.vector("community_rate", std::move(commu));
        out.matrix("p2p_rate", std::move(p2p));
    }

    //////////////////////////////  Model : transmission cache  //////////////////////////////

    /// Store profiles between infected members `indiv_i` and `indiv_j`, by generation time.
    void helper_update_pairwise_cache(matrix &inst_mtx,
                                      matrix &cumul_mtx,
                                      const Individual &indiv_i,
                                      const Individual &indiv_j,
                                      double t_infec_i,
                                      double t_infec_j)
    {
        const MemberIdx i = Index::member(indiv_i);
        const MemberIdx j = Index::member(indiv_j);
        const double gentime = t_infec_j - t_infec_i;

        // the reverse direction is set to 0: the rate of infection sums over every pair
        if (gentime > 0)
        { // i : infector, j : infectee
            inst_mtx[i][j] = infec_profile_inst(indiv_i, gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
            cumul_mtx[i][j] = infec_profile_cumul(indiv_i, gentime - 1);
#else
            cumul_mtx[i][j] = infec_profile_cumul(indiv_i, gentime);
#endif
            inst_mtx[j][i] = 0.0;
            cumul_mtx[j][i] = 0.0;
        }
        else if (gentime < 0)
        { // j : infector, i : infectee
            inst_mtx[i][j] = 0.0;
            cumul_mtx[i][j] = 0.0;
            inst_mtx[j][i] = infec_profile_inst(indiv_j, -gentime);
#ifdef _HOUSEPI_DISCRETE_TIME
            cumul_mtx[j][i] = infec_profile_cumul(indiv_j, -gentime - 1);
#else
            cumul_mtx[j][i] = infec_profile_cumul(indiv_j, -gentime);
#endif
        }
#ifdef _HOUSEPI_DISCRETE_TIME
        else
        {
            inst_mtx[i][j] = 0.0;
            cumul_mtx[i][j] = 0.0;
            inst_mtx[j][i] = 0.0;
            cumul_mtx[j][i] = 0.0;
        }
#endif
    }

    /// Commit the introduction: earliest infection, ties to the smallest member index.
    static void helper_commit_intro(MemberIdx &i0,
                                    double &t0,
                                    MemberIdx i,
                                    double t)
    {
        if (t < t0 || (t == t0 && i < i0))
        {
            t0 = t;
            i0 = i;
        }
    }

    /// Store cumulative and instantaneous infectivity profiles of `house`.
    /// @note Must not read a mutable parameter.
    void update_house_transmission_cache(const Household &house) override final
    {
        auto &inst_mtx = cache_double_matrix(house, HouseDoubleMtx::inst);
        auto &cumul_mtx = cache_double_matrix(house, HouseDoubleMtx::cumul);
        MemberIdx &i0 = cache_member(house, HouseMember::intro_case);
        double new_t0 = 100000;

        const MemberIdx n_infec = house.n_infected_status();
        const auto &infected_members = house.infected_members();

        for (MemberIdx a = 0; a < n_infec; ++a)
        {
            const MemberIdx i = infected_members[a];
            const Individual &indiv_i = house.at(i);
            const double t_infec_i = aug(Latent::infec_time, indiv_i);

            //////  Escaped members  //////

            for (MemberIdx s : house.escaped_members())
            {
                const Individual &indiv_s = house.at(s);
                const double exposure = follow_up_end(indiv_s) - t_infec_i;
                cumul_mtx[i][s] = infec_profile_cumul(indiv_i, exposure);
                cumul_mtx[s][i] = 0.0;
            }

            //////  Infectious members  //////

            for (MemberIdx b = a + 1; b < n_infec; ++b)
            {
                const MemberIdx j = infected_members[b];
                const Individual &indiv_j = house.at(j);
                const double t_infec_j = aug(Latent::infec_time, indiv_j);

                helper_update_pairwise_cache(inst_mtx, cumul_mtx,
                                             indiv_i, indiv_j,
                                             t_infec_i, t_infec_j);
            }

            // update household introduction case
            helper_commit_intro(i0, new_t0, i, t_infec_i);
        }
    }

    /// Store the infectivity profiles involving infector `indiv_i`.
    /// @note Undefined when `indiv_i` is not an infector.
    void update_member_transmission_cache(const Household &house,
                                          const Individual &indiv_i) override final
    {
        auto &inst_mtx = cache_double_matrix(house, HouseDoubleMtx::inst);
        auto &cumul_mtx = cache_double_matrix(house, HouseDoubleMtx::cumul);
        MemberIdx &i0 = cache_member(house, HouseMember::intro_case);

        const MemberIdx i = Index::member(indiv_i);
        const double t_infec_i = aug(Latent::infec_time, indiv_i);

        //////  Escaped members  //////

        for (MemberIdx s : house.escaped_members())
        {
            const Individual &indiv_s = house.at(s);
            double exposure = follow_up_end(indiv_s) - t_infec_i;
            cumul_mtx[i][s] = infec_profile_cumul(indiv_i, exposure);
            cumul_mtx[s][i] = 0.0;
        }

        //////  Infectious members  //////

        // commit introduction to current member
        i0 = i;
        double new_t0 = t_infec_i;

        for (MemberIdx j : house.infected_members())
        {
            if (j == i)
                continue;

            const Individual &indiv_j = house.at(j);
            const double t_infec_j = aug(Latent::infec_time, indiv_j);

            helper_update_pairwise_cache(inst_mtx, cumul_mtx,
                                         indiv_i, indiv_j,
                                         t_infec_i, t_infec_j);

            // update household introduction case
            helper_commit_intro(i0, new_t0, j, t_infec_j);
        }
    }

    //////////////////////////////  Model : household simulation  //////////////////////////////

    /// Simulate the outbreak in `house` as a Poisson point process.
    void simulate_single_house(Household &house) override final
    {
#ifdef _HOUSEPI_DISCRETE_TIME
        const double Dt = 1;
#else
        const double Dt = real_time_step();
#endif
        //////  Generation 0: single introduction  //////

        const MemberIdx i0 = set_introduction(house, Latent::infec_time);
        Individual &indiv_0 = house.at(i0);

        //////  Generation 1+: infection by the previous generation  //////

        const double alpha = par(Param::alpha);
        const double beta = transmission_rate(house.size());

        const double t0 = aug(Latent::infec_time, indiv_0);
        for (double t = t0 + Dt; house.n_escaped_status() > 0; t += Dt)
        {
            bool any_at_risk = false;
            for (Individual &indiv_i : house)
            {
                if (indiv_i.is_infected_status())
                    continue;
                if (t > follow_up_end(indiv_i))
                    continue;
                any_at_risk = true;

                if (exponential() < roi_inst_simulation(house, indiv_i, t, alpha, beta) * Dt)
                {
                    set_infected(house, indiv_i);
                    set_latent_value(Latent::infec_time, indiv_i, t);
                    init_missing_latent(indiv_i);
                }
            }
            if (!any_at_risk)
                break;
        }
    }

    /// Instantaneous rate of infection of `indiv_i` at time `t_infec_i`, in simulation.
    double roi_inst_simulation(
        const Household &house,
        const Individual &indiv_i,
        double t_infec_i,
        double alpha,
        double beta) const
    {
        MemberIdx i = Index::member(indiv_i);
        double roi = 0.; // instantaneous rate of infection

        // j : infector, i : infectee
        for (MemberIdx j : house.infected_members())
        {
            const Individual &indiv_j = house.at(j);
            const double gentime = t_infec_i - aug(Latent::infec_time, indiv_j);
            roi += group_rate(GroupRate::infectivity, indiv_j) *
                   contact_rate(house, j, i) *
                   infec_profile_inst(indiv_j, gentime);
        }

        return alpha * group_rate(GroupRate::community, indiv_i) +
               beta * group_rate(GroupRate::susceptibility, indiv_i) * roi;
    }
};
HOUSEPI_CLOSE_MODULE_NAMESPACE()
