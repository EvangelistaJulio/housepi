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
	HOUSEPI_DECLARE_EMPTY_KEYS(HHsizeKey)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHsizevecKey, order, types)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublevecKey, commu_escape_rate)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublemtxKey, p2p_escape_rate)
	HOUSEPI_DECLARE_REGISTER_KEYS(ParamKey, alpha, beta)

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

		//////////////////////////////////
		/////////////  TYPES  ////////////
		//////////////////////////////////

		// initialize types for house escape matrix
		// NOTE: cache vector initialized to 0
		auto &types = cache_size_vector(house, HHsizevecKey::types);

		// identify distinct types among seropositives
		// by checking known group rate attributes direct comparison
		// UPDATE: can't I use DataWide multiplicituy mechanism instead?
		size_t n_infec = house.n_infected_status();
		std::vector<const Attribute *> community_attr(n_infec);
		std::vector<const Attribute *> infectivity_attr(n_infec);
		std::vector<const Attribute *> susceptibility_attr(n_infec);

		// progressive: first infected type always 1, then 2, etc. seronegative always type 0
		size_t type_counter = 0;
		for (size_t i = 0; i < n_infec; ++i)
		{
			const Individual &member_i = house.member(order[i]);
			community_attr[i] = &member_i.get_group_rate_attribute(GroupRateKey::community);
			infectivity_attr[i] = &member_i.get_group_rate_attribute(GroupRateKey::infectivity);
			susceptibility_attr[i] = &member_i.get_group_rate_attribute(GroupRateKey::susceptibility);

			++type_counter;
			types[i] = type_counter;
			for (size_t j = 0; j < i; ++j)
				if (*community_attr[j] == *community_attr[i] &&
					*infectivity_attr[j] == *infectivity_attr[i] &&
					*susceptibility_attr[j] == *susceptibility_attr[i])
				{
					--type_counter;
					types[i] = types[j];
					break;
				}
		}
	}

	/// Risk of infection from the community.
	double community_risk() const { return par(ParamKey::alpha); }
	/// Risk of infection from another household member.
	double person_to_person_risk() const { return par(ParamKey::beta); }
	/// Store person-to-person escape probabilities inside a cache.
	void fill_house_escape_matrix(Household &house);

	/// Compute the log-likelihood for a single household.
	double compute_house_obs_log_lik(Household &house) override final
	{
		// fill per-individual escape rate matrices
		fill_house_escape_matrix(house);

		return bino::chain_binom_grouped_log_lik(
			house.n_infected_status(),
			cache_size_vector(house, HHsizevecKey::types),
			cache_double_vector(house, HHdoublevecKey::commu_escape_rate),
			cache_double_matrix(house, HHdoublemtxKey::p2p_escape_rate));
	}

	/// Simulate the infection within a household.
	void simulate_single_house(Household &house) override final;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Model : compute single house loglike  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::fill_house_escape_matrix(Household &house)
{
	const auto &order = cache_size_vector(house, HHsizevecKey::order);
	const double commu_rate = community_risk();
	const double p2p_rate = person_to_person_risk();

	// assign values to escape matrix
	auto &commu_mtx = cache_double_vector(house, HHdoublevecKey::commu_escape_rate);
	auto &p2p_mtx = cache_double_matrix(house, HHdoublemtxKey::p2p_escape_rate);

	// fill matrix according to order
	for (size_t i = 0; i < house.size(); ++i)
	{
		const Individual &member_i = house.member(order[i]);
		double inf_i = member_i.group_rate(GroupRateKey::infectivity);
		double sus_i = member_i.group_rate(GroupRateKey::susceptibility);

		// community
		commu_mtx[i] = commu_rate * member_i.group_rate(GroupRateKey::community);

		// person-to-person
		// NOTE: with respect to heterogeneous, I need to compute also diagonal terms
		// UPDATE: I would like a unifying approach...
		for (size_t k = 0; k <= i; ++k)
		{
			const Individual &member_k = house.member(order[k]);
			// [i] : infector -> [k] : infectee
			p2p_mtx[i][k] = inf_i * p2p_rate * member_k.group_rate(GroupRateKey::susceptibility);
			// [k] : infector -> [i] : infectee
			p2p_mtx[k][i] = member_k.group_rate(GroupRateKey::infectivity) * p2p_rate * sus_i;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::simulate_single_house(Household &house)
{
	const double commu_rate = community_risk();
	const double p2p_rate = person_to_person_risk();

	//////////////////////////////////////////////////////////////////
	/////////////  Generation 0 : infection by community  ////////////
	//////////////////////////////////////////////////////////////////

	std::vector<Individual *> new_gen_set = {};

	for (Individual &member : house)
		if (exponential() > commu_rate * member.group_rate(GroupRateKey::community))
			house.set_susceptible_status(member);
		else
		{
			house.set_infected_status(member);
			new_gen_set.push_back(&member);
		}

	//////////////////////////////////////////////////////////////////
	//////  Generation >=1 : infection by previous generation  ///////
	//////////////////////////////////////////////////////////////////

	std::vector<Individual *> old_gen_set = {};

	while (house.n_susceptible_status() > 0 && !new_gen_set.empty())
	{
		std::swap(old_gen_set, new_gen_set);
		new_gen_set.clear();

		for (Individual &member_i : house)
		{
			if (member_i.is_infected_status())
				continue;

			double r = 0.;
			for (Individual *member_k : old_gen_set)
				r += member_k->group_rate(GroupRateKey::infectivity);

			double sus_i = member_i.group_rate(GroupRateKey::susceptibility);
			if (exponential() < r * p2p_rate * sus_i)
			{
				house.set_infected_status(member_i);
				new_gen_set.push_back(&member_i);
			}
		}
	}
}
HOUSEPI_CLOSE_MODULE_NAMESPACE()
