/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * CHAIN BINOMIAL * *** *  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "housepi.h"
#include "binomial.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE final_size
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()
class Model final : public ModelBase
{
public:
	HOUSEPI_DECLARE_EMPTY_KEYS(Latent)
	HOUSEPI_DECLARE_EMPTY_KEYS(Data)
	HOUSEPI_DECLARE_EMPTY_KEYS(Long)
	HOUSEPI_DECLARE_REGISTER_KEYS(GroupRisk, community, susceptibility, infectivity)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHmember, num_parts)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHmembervec, order, part_counts, part_indexes)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublevec, commu_escape_rate)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublemtx, p2p_escape_rate)
	HOUSEPI_DECLARE_REGISTER_KEYS(Param, alpha, beta, gamma, hhsize_ref)

	/// General advanced inscrition for the setup.
	void init_model_requirements() override final
	{
		request_param_discrete(
			{Param::hhsize_ref});
		request_risk_identifiability(
			Param::alpha, {GroupRisk::community});
		request_risk_identifiability(
			Param::beta, {GroupRisk::susceptibility, GroupRisk::infectivity}, true);
	}

	/// Allow for contact risk definitions.
	bool has_contact_risk() const override final { return true; }

	void init_model_for_evaluation() override final
	{
		assert_no_missing_status();
		chain_binom::clear_lookuptables();
	}

	/// Initialize hosuehold for evaluation
	void init_house_for_evaluation(const Household &house) override final
	{
		//////////////////////////////////
		/////////////  ORDER  ////////////
		//////////////////////////////////

		// initialize internal order for house escape matrix
		// ordered members : 1) infected, 2) escaped
		auto &order = cache_member_vector(house, HHmembervec::order);
		order.clear();
		order.reserve(house.size());
		order.insert(order.end(), house.infected_members().begin(), house.infected_members().end());
		order.insert(order.end(), house.escaped_members().begin(), house.escaped_members().end());

		///////////////////////////////////////
		/////////////  PARTITIONS  ////////////
		///////////////////////////////////////

		// associate infected integer-partition, built in place in the cache:
		// both vectors hold house.size() entries and K never exceeds n_infec
		const MemberIndex n_infec = house.n_infected_status();
		auto &part_counts = cache_member_vector(house, HHmembervec::part_counts);
		auto &part_indexes = cache_member_vector(house, HHmembervec::part_indexes);
		MemberIndex K = 0; // number of partitions

		for (MemberIndex oi = 0; oi < n_infec; ++oi)
		{
			// search existing partitions for a match
			MemberIndex p = K; // tentatively new
			for (MemberIndex oj = 0; oj < K; ++oj)
				if (is_same_risk_partition(house, order[part_indexes[oj]], order[oi]))
					p = oj; // this acts as a break

			if (p == K)
			{
				// append new partition at the end
				part_counts[K] = 1;
				part_indexes[K] = oi;
				++K;
			}
			else
			{
				// update counts, then bubble-up to guarantee descending order
				++part_counts[p];
				while (p > 0 && part_counts[p] > part_counts[p - 1])
				{
					std::swap(part_counts[p], part_counts[p - 1]);
					std::swap(part_indexes[p], part_indexes[p - 1]);
					--p;
				}
			}
		}

		cache_member(house, HHmember::num_parts) = K;

		if (K <= 1 && n_infec > 1)
			chain_binom::ensure_homogeneous(n_infec);
		else if (K > 1 && K < n_infec)
			chain_binom::ensure_partitioned(house.eq_class(), part_counts, K);
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////  Model : compute single house loglike  ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	/// Risk of infection from the community.
	double community_risk() const { return par(Param::alpha); }
	/// Risk of infection from another household member.
	double person_to_person_risk(MemberIndex house_size) const
	{
		return par(Param::beta) *
			   std::exp(par(Param::gamma) * std::log(par(Param::hhsize_ref) / house_size));
	}

	/// Store person-to-person escape probabilities inside a cache.
	/// UPDATE: I don't need to compute for all members, just for partitions.
	/// This requires to compute also escaped partitions.
	void fill_house_escape_matrix(const Household &house)
	{
		const MemberIndex house_size = house.size();
		const auto &order = cache_member_vector(house, HHmembervec::order);
		const double commu_rate = community_risk();
		const double p2p_rate = person_to_person_risk(house_size);

		// assign values to escape matrix
		auto &commu_mtx = cache_double_vector(house, HHdoublevec::commu_escape_rate);
		auto &p2p_mtx = cache_double_matrix(house, HHdoublemtx::p2p_escape_rate);

		// fill matrix according to order
		for (MemberIndex oi = 0; oi < house_size; ++oi)
		{
			MemberIndex i = order[oi]; // retrieve household member index (partition representative)
			const Individual &indiv_i = house.at(i);
			const double inf_i_scl = group_risk(GroupRisk::infectivity, indiv_i) * p2p_rate;
			const double sus_i = group_risk(GroupRisk::susceptibility, indiv_i);
			const double sus_i_scl = p2p_rate * sus_i;

			// community filled in order
			commu_mtx[oi] = commu_rate * group_risk(GroupRisk::community, indiv_i);

			// person-to-person
			for (MemberIndex oj = 0; oj < oi; ++oj)
			{
				MemberIndex j = order[oj];
				const Individual &indiv_j = house.at(j);
				// [i] : infector -> [k] : infectee
				p2p_mtx[oi][oj] = inf_i_scl * contact_risk(house, i, j) *
								  group_risk(GroupRisk::susceptibility, indiv_j);
				// [k] : infector -> [i] : infectee
				p2p_mtx[oj][oi] = group_risk(GroupRisk::infectivity, indiv_j) *
								  contact_risk(house, j, i) * sus_i_scl;
			}

			// NOTE: for partitions, I need to compute also diagonal terms:
			// i. e. from partition 'i' to partition 'i' (intra-partition)
			// NOTE: by definition contact_risk(house, i, i) = 0
			p2p_mtx[oi][oi] = inf_i_scl * intra_partition_contact_risk(house, i) * sus_i;
		}
	}

	/// Compute transmission log-likelihood within household.
	double compute_house_transm_log_lik(const Household &house) override final
	{
		// fill per-individual escape rate matrices
		fill_house_escape_matrix(house);

		return chain_binom::log_likelihood(
			house.eq_class(),
			house.n_infected_status(),
			cache_member(house, HHmember::num_parts),
			cache_member_vector(house, HHmembervec::part_indexes),
			cache_double_vector(house, HHdoublevec::commu_escape_rate),
			cache_double_matrix(house, HHdoublemtx::p2p_escape_rate));
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////  Model : simulate single house  //////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	/// Simulate the infection within a household.
	void simulate_single_house(Household &house) override final
	{
		const double commu_rate = community_risk();
		const double p2p_rate = person_to_person_risk(house.size());

		//////////////////////////////////////////////////////////////////
		/////////////  Generation 0 : infection by community  ////////////
		//////////////////////////////////////////////////////////////////

		std::vector<MemberIndex> new_gen_set = {};

		for (Individual &indiv : house)
			if (exponential() > commu_rate * group_risk(GroupRisk::community, indiv))
				set_escaped(house, indiv);
			else
			{
				set_infected(house, indiv);
				new_gen_set.push_back(Index::member(indiv));
			}

		//////////////////////////////////////////////////////////////////
		//////  Generation >=1 : infection by previous generation  ///////
		//////////////////////////////////////////////////////////////////

		std::vector<MemberIndex> old_gen_set = {};

		while (house.n_escaped_status() > 0 && !new_gen_set.empty())
		{
			std::swap(old_gen_set, new_gen_set);
			new_gen_set.clear();

			for (MemberIndex i = 0; i < house.size(); ++i)
			{
				Individual &indiv_i = house.at(i);

				if (indiv_i.is_infected_status())
					continue;

				double r = 0.;
				for (MemberIndex j : old_gen_set)
					r += group_risk(GroupRisk::infectivity, house.at(j)) *
						 contact_risk(house, j, i);

				double sus_i = group_risk(GroupRisk::susceptibility, indiv_i);
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
