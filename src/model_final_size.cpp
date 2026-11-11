/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

//////////////////////////////  final_size  //////////////////////////////

#include <memory>

#include "chain_binom.h"
#include "housepi.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE final_size
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()
/// Model final-size outbreaks with a chain binomial: community and person-to-person escape.
class Model final : public ModelBase
{
private:
    /// Chain-binomial tables, shared by the copies of the model.
    std::shared_ptr<chain_binom::LookUpTables> lut_ = std::make_shared<chain_binom::LookUpTables>();
    /// Shape index in `lut_`, per class representative household.
    std::vector<chain_binom::ShapeIdx> shape_idx_;
    /// Chain-binomial buffers, per thread slot.
    std::vector<chain_binom::Scratch> scratch_;

public:
    HOUSEPI_DECLARE_GROUP_RATE_KEYS(community, susceptibility, infectivity)
    HOUSEPI_DECLARE_HOUSE_MEMBER_KEYS(num_parts)
    HOUSEPI_DECLARE_HOUSE_MEMBER_VEC_KEYS(order, part_counts, part_indexes)
    HOUSEPI_DECLARE_HOUSE_DOUBLE_VEC_KEYS(commu_escape_rate)
    HOUSEPI_DECLARE_HOUSE_DOUBLE_MTX_KEYS(p2p_escape_rate)
    HOUSEPI_DECLARE_PARAM_KEYS(alpha, beta, gamma, hh_size_ref)

    /// Declare parameter and rate requirements.
    void init_model_requirements() override final
    {
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
    }

    /// Whether the model uses contact rates.
    bool has_contact_rate() const override final { return true; }

    /// Assert no member has missing status, then allocate buffers per thread slot.
    void init_model_for_evaluation() override final
    {
        assert_no_missing_status();
        scratch_.resize(n_thread_slots());
    }

    /// Order members of `house`, then group infected ones into rate partitions.
    /// * partitions needing more than chain_binom::MAX_TERMS terms -> error.
    void init_house_for_evaluation(const Household &house) override final
    {
        //////  Order  //////

        // order of the escape matrix: infected members first, then escaped ones
        auto &order = cache_member_vector(house, HouseMemberVec::order);
        order.clear();
        order.reserve(house.size());
        order.insert(order.end(), house.infected_members().begin(), house.infected_members().end());
        order.insert(order.end(), house.escaped_members().begin(), house.escaped_members().end());

        //////  Partitions  //////

        // integer partition of infected members, built in place in the cache
        // both vectors hold house.size() entries, K never exceeds n_infec
        const MemberIdx n_infec = house.n_infected_status();
        auto &part_counts = cache_member_vector(house, HouseMemberVec::part_counts);
        auto &part_indexes = cache_member_vector(house, HouseMemberVec::part_indexes);
        // number of partitions
        MemberIdx K = 0;

        for (MemberIdx oi = 0; oi < n_infec; ++oi)
        {
            // search existing partitions for a match
            // tentatively new
            MemberIdx p = K;
            for (MemberIdx oj = 0; oj < K; ++oj)
                if (is_same_rate_partition(house, order[part_indexes[oj]], order[oi]))
                    // partitions are disjoint: at most one matches
                    p = oj;

            if (p == K)
            {
                // append new partition at the end
                part_counts[K] = 1;
                part_indexes[K] = oi;
                ++K;
            }
            else
            {
                // update count, then bubble up to keep counts descending
                ++part_counts[p];
                while (p > 0 && part_counts[p] > part_counts[p - 1])
                {
                    std::swap(part_counts[p], part_counts[p - 1]);
                    std::swap(part_indexes[p], part_indexes[p - 1]);
                    --p;
                }
            }
        }
        cache_member(house, HouseMember::num_parts) = K;

        // assert that the model is tractable
        chain_binom::assert_tractable(part_counts, K, str(house));

        // assign a shape index to the class representative
        if (shape_idx_.size() <= Index::of(house))
            shape_idx_.resize(Index::of(house) + 1);
        shape_idx_[Index::of(house)] = (K > 1 && K < n_infec) ? lut_->ensure_partitioned(part_counts, K) : 0;
    }

    //////////////////////////////  Model : household log-likelihood  //////////////////////////////

    /// Return rate of infection from the community.
    double community_rate() const { return par(Param::alpha); }
    /// Return rate of infection from another member, scaled by household size.
    /// Based on @cite cauchemez_bayesian_2004
    double person_to_person_rate(MemberIdx house_size) const
    {
        return par(Param::beta) *
               std::exp(par(Param::gamma) * std::log(par(Param::hh_size_ref) / house_size));
    }

    /// Store community and person-to-person escape rates of `house` in the cache.
    /// @todo Compute per partition, not per member, escaped partitions included:
    /// * fill rates between part representatives (`part_indexes`) only;
    /// * partition escaped members in `init_house_for_evaluation()`;
    /// * sum escaped contributions by counts, sum_E n_E (commu_E + sum_P n_P p2p_PE);
    /// * read the intra-partition diagonal in the homogeneous branch, not `[0][1]`;
    /// * keep a per-member fill for export_house_quantities().
    void fill_house_escape_matrix(const Household &house)
    {
        const MemberIdx house_size = house.size();
        const auto &order = cache_member_vector(house, HouseMemberVec::order);
        const double commu_rate = community_rate();
        const double p2p_rate = person_to_person_rate(house_size);

        // escape rates
        auto &commu_mtx = cache_double_vector(house, HouseDoubleVec::commu_escape_rate);
        auto &p2p_mtx = cache_double_matrix(house, HouseDoubleMtx::p2p_escape_rate);

        // fill in member order
        for (MemberIdx oi = 0; oi < house_size; ++oi)
        {
            // member index, partition representative
            MemberIdx i = order[oi];
            const Individual &indiv_i = house.at(i);
            const double inf_i_scl = group_rate(GroupRate::infectivity, indiv_i) * p2p_rate;
            const double sus_i = group_rate(GroupRate::susceptibility, indiv_i);
            const double sus_i_scl = p2p_rate * sus_i;

            // community
            commu_mtx[oi] = commu_rate * group_rate(GroupRate::community, indiv_i);

            // person-to-person
            for (MemberIdx oj = 0; oj < oi; ++oj)
            {
                MemberIdx j = order[oj];
                const Individual &indiv_j = house.at(j);
                // i infector, j infectee
                p2p_mtx[oi][oj] = inf_i_scl * contact_rate(house, i, j) *
                                  group_rate(GroupRate::susceptibility, indiv_j);
                // j infector, i infectee
                p2p_mtx[oj][oi] = group_rate(GroupRate::infectivity, indiv_j) *
                                  contact_rate(house, j, i) * sus_i_scl;
            }

            // diagonal: intra-partition term, read by partitions
            // contact_rate(house, i, i) is 0 by definition
            p2p_mtx[oi][oi] = inf_i_scl * intra_partition_contact_rate(house, i) * sus_i;
        }
    }

    /// Compute transmission log-likelihood within household.
    double compute_house_transm_log_lik(const Household &house) override final
    {
        // per-member escape rates
        fill_house_escape_matrix(house);

        return chain_binom::log_likelihood(
            *lut_,
            scratch_[thread_slot()],
            house.n_infected_status(),
            cache_member(house, HouseMember::num_parts),
            shape_idx_[Index::of(house)],
            cache_member_vector(house, HouseMemberVec::part_indexes),
            cache_double_vector(house, HouseDoubleVec::commu_escape_rate),
            cache_double_matrix(house, HouseDoubleMtx::p2p_escape_rate));
    }

    /// Export community and person-to-person escape probabilities of `house`, by member.
    /// * `community_escape`: one per member;
    /// * `p2p_escape`: infector by row, infectee by column, diagonal 1.
    void export_house_quantities(const Household &house, HouseQuantities &out) override final
    {
        // escape rates at the current parameters
        fill_house_escape_matrix(house);

        const MemberIdx house_size = house.size();
        const auto &order = cache_member_vector(house, HouseMemberVec::order);
        const auto &commu_mtx = cache_double_vector(house, HouseDoubleVec::commu_escape_rate);
        const auto &p2p_mtx = cache_double_matrix(house, HouseDoubleMtx::p2p_escape_rate);

        // from member order back to member index
        std::vector<double> commu(house_size);
        matrix p2p(house_size, std::vector<double>(house_size, 1.));
        for (MemberIdx oi = 0; oi < house_size; ++oi)
        {
            commu[order[oi]] = std::exp(-commu_mtx[oi]);
            for (MemberIdx oj = 0; oj < house_size; ++oj)
                if (oj != oi)
                    p2p[order[oi]][order[oj]] = std::exp(-p2p_mtx[oi][oj]);
        }
        out.vector("community_escape", std::move(commu));
        out.matrix("p2p_escape", std::move(p2p));
    }

    //////////////////////////////  Model : household simulation  //////////////////////////////

    /// Simulate transmission within `house`, generation by generation.
    void simulate_single_house(Household &house) override final
    {
        const double commu_rate = community_rate();
        const double p2p_rate = person_to_person_rate(house.size());

        //////  Generation 0: community infection  //////

        std::vector<MemberIdx> new_gen_set = {};

        for (Individual &indiv : house)
            if (exponential() > commu_rate * group_rate(GroupRate::community, indiv))
                set_escaped(house, indiv);
            else
            {
                set_infected(house, indiv);
                new_gen_set.push_back(Index::member(indiv));
            }

        //////  Generation 1+: infection by the previous generation  //////

        std::vector<MemberIdx> old_gen_set = {};

        while (house.n_escaped_status() > 0 && !new_gen_set.empty())
        {
            std::swap(old_gen_set, new_gen_set);
            new_gen_set.clear();

            for (MemberIdx i = 0; i < house.size(); ++i)
            {
                Individual &indiv_i = house.at(i);

                if (indiv_i.is_infected_status())
                    continue;

                double r = 0.;
                for (MemberIdx j : old_gen_set)
                    r += group_rate(GroupRate::infectivity, house.at(j)) *
                         contact_rate(house, j, i);

                double sus_i = group_rate(GroupRate::susceptibility, indiv_i);
                if (exponential() < r * p2p_rate * sus_i)
                {
                    set_infected(house, indiv_i);
                    new_gen_set.push_back(i);
                }
            }
        }
    }
};
HOUSEPI_CLOSE_MODULE_NAMESPACE()
