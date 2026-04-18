/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  * *** * CHAIN BINOMIAL GROUPED * *** *  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "housepi.h"
#include "binomial.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE chain_binomial_grouped
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()
class Model final : public ModelBase
{
public:
	HOUSEPI_DECLARE_EMPTY_KEYS(AugmKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
	HOUSEPI_DECLARE_REGISTER_KEYS(GroupRateKey, community, susceptibility, infectivity)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHsizeKey, num_parts)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHsizevecKey, order, part_counts, part_indexes)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublevecKey, commu_escape_rate)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublemtxKey, p2p_escape_rate)
	HOUSEPI_DECLARE_REGISTER_KEYS(ParamKey, alpha, beta, gamma, size_ref)

	/// General advanced inscrition for the setup.
	void init_advanced_setup() override final
	{
		request_discrete_param(ParamKey::size_ref);
	}

	/// Allow for contact rate definitions.
	bool has_contact_rates() const override final { return true; }

	/// Initialize the LookUp Tables to largest possible values.
	/// NOTE: missing status are not modeled.
	void init_model_for_evaluation() override final
	{
		assert_no_missing_status();
		/// UPDATE: is this actually used?
		bino::init_generations_lookuptable(get_max_hh_cases());
	}

	/// Initialize hosuehold for evaluation
	void init_house_for_evaluation(Household &house) override final
	{
		//////////////////////////////////
		/////////////  ORDER  ////////////
		//////////////////////////////////

		// initialize internal order for house escape matrix
		// ordered members : 1) infectious, 2) susceptibles
		auto &order = cache_size_vector(house, HHsizevecKey::order);
		order.clear();
		order.reserve(house.size());
		order.insert(order.end(), house.infected_hhindexes().begin(), house.infected_hhindexes().end());
		order.insert(order.end(), house.susceptible_hhindexes().begin(), house.susceptible_hhindexes().end());

		///////////////////////////////////////
		/////////////  PARTITIONS  ////////////
		///////////////////////////////////////

		// associate infected integer-partition
		const size_t n_infec = house.n_infected_status();
		std::vector<size_t> partition_indexes = {};	 // representatives
		std::vector<unsigned> partition_counts = {}; // counts

		// UPDATE: I will move this in a ModelBase helper function
		for (size_t oi = 0; oi < n_infec; ++oi)
		{
			// search existing partitions for a match
			size_t K = partition_indexes.size(); // tentatively new
			for (size_t oj = 0; oj < K; ++oj)
				if (is_same_rate_partition(house, order[partition_indexes[oj]], order[oi]))
					K = oj; // this acts as a break

			if (K == partition_indexes.size())
			{
				// append new partition at the end
				partition_counts.push_back(1);
				partition_indexes.push_back(oi);
			}
			else
			{
				// update counts, then bubble-up to guarantee descending order
				++partition_counts[K];
				while (K > 0 && partition_counts[K] > partition_counts[K - 1])
				{
					std::swap(partition_counts[K], partition_counts[K - 1]);
					std::swap(partition_indexes[K], partition_indexes[K - 1]);
					--K;
				}
			}
		}

		// store in the cache
		size_t K = partition_indexes.size(); // n. of infected partitions
		for (size_t i = 0; i < K; ++i)
		{
			cache_size_vector(house, HHsizevecKey::part_counts)[i] = partition_counts[i];
			cache_size_vector(house, HHsizevecKey::part_indexes)[i] = partition_indexes[i];
		}
		cache_size(house, HHsizeKey::num_parts) = K;
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////  Model : compute single house loglike  ////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////////

	/// Risk of infection from the community.
	double community_risk() const { return par(ParamKey::alpha); }
	/// Risk of infection from another household member.
	double person_to_person_risk(size_t house_size) const
	{
		return par(ParamKey::beta) *
			   std::exp(par(ParamKey::gamma) * std::log(par(ParamKey::size_ref) / house_size));
	}

	/// Store person-to-person escape probabilities inside a cache.
	/// UPDATE: I don't need to compute for all members, just for partitions.
	/// This requires to compute also susceptible partitions.
	void fill_house_escape_matrix(Household &house)
	{
		size_t house_size = house.size();
		const auto &order = cache_size_vector(house, HHsizevecKey::order);
		const double commu_rate = community_risk();
		const double p2p_rate = person_to_person_risk(house_size);

		// assign values to escape matrix
		auto &commu_mtx = cache_double_vector(house, HHdoublevecKey::commu_escape_rate);
		auto &p2p_mtx = cache_double_matrix(house, HHdoublemtxKey::p2p_escape_rate);

		// fill matrix according to order
		for (size_t oi = 0; oi < house_size; ++oi)
		{
			size_t i = order[oi]; // retrieve household member index (partition representative)
			const Individual &member_i = house.member(i);
			const double inf_i_scl = member_i.group_rate(GroupRateKey::infectivity) * p2p_rate;
			const double sus_i = member_i.group_rate(GroupRateKey::susceptibility);
			const double sus_i_scl = p2p_rate * sus_i;

			// community filled in order
			commu_mtx[oi] = commu_rate * member_i.group_rate(GroupRateKey::community);

			// person-to-person
			for (size_t oj = 0; oj < oi; ++oj)
			{
				size_t j = order[oj];
				const Individual &member_j = house.member(j);
				// [i] : infector -> [k] : infectee
				p2p_mtx[oi][oj] = inf_i_scl * house.contact_rate(i, j) *
								  member_j.group_rate(GroupRateKey::susceptibility);
				// [k] : infector -> [i] : infectee
				p2p_mtx[oj][oi] = member_j.group_rate(GroupRateKey::infectivity) *
								  house.contact_rate(j, i) * sus_i_scl;
			}

			// NOTE: four partitions, I need to compute also diagonal terms:
			// i. e. from partition 'i' to partition 'i' (intra-partition)
			// NOTE: by definition house.contact_rate(i, i) = 0
			p2p_mtx[oi][oi] = inf_i_scl * house.intra_partition_contact_rate(i) * sus_i;
		}
	}

	/// Compute the log-likelihood for a single household.
	double compute_house_obs_log_lik(Household &house) override final
	{
		// fill per-individual escape rate matrices
		fill_house_escape_matrix(house);

		return bino::chain_binom_grouped_log_lik(
			house.n_infected_status(),
			cache_size(house, HHsizeKey::num_parts),
			cache_size_vector(house, HHsizevecKey::part_counts),
			cache_size_vector(house, HHsizevecKey::part_indexes),
			cache_double_vector(house, HHdoublevecKey::commu_escape_rate),
			cache_double_matrix(house, HHdoublemtxKey::p2p_escape_rate));
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

		std::vector<size_t> new_gen_set = {};

		for (size_t i = 0; i < house.size(); ++i)
			if (exponential() > commu_rate * house.member(i).group_rate(GroupRateKey::community))
				house.set_susceptible_status(i);
			else
			{
				house.set_infected_status(i);
				new_gen_set.push_back(i);
			}

		//////////////////////////////////////////////////////////////////
		//////  Generation >=1 : infection by previous generation  ///////
		//////////////////////////////////////////////////////////////////

		std::vector<size_t> old_gen_set = {};

		while (house.n_susceptible_status() > 0 && !new_gen_set.empty())
		{
			std::swap(old_gen_set, new_gen_set);
			new_gen_set.clear();

			for (size_t i = 0; i < house.size(); ++i)
			{
				Individual &member_i = house.member(i);

				if (member_i.is_infected_status())
					continue;

				double r = 0.;
				for (size_t j : old_gen_set)
					r += house.member(j).group_rate(GroupRateKey::infectivity) *
						 house.contact_rate(j, i);

				double sus_i = member_i.group_rate(GroupRateKey::susceptibility);
				if (exponential() < r * p2p_rate * sus_i)
				{
					house.set_infected_status(i);
					new_gen_set.push_back(i);
				}
			}
		}
	}
};
HOUSEPI_CLOSE_MODULE_NAMESPACE()
