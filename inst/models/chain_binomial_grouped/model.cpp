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

#include <algorithm>
#include <iterator>
#include <map>
#include <tuple>

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
	HOUSEPI_DECLARE_REGISTER_KEYS(HHsizevecKey, order)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublevecKey, commu_escape_rate)
	HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublemtxKey, p2p_escape_rate)
	HOUSEPI_DECLARE_REGISTER_KEYS(ParamKey, alpha, beta)

	/// Initialize the LookUp Tables to largest possible values.
	/// NOTE: missing status are not modeled.
	void init_model_for_evaluation() override final
	{
		assert_no_missing_status();
		bino::init_generations_lookuptable(get_max_hh_cases());
	}

	/// Initialize internal order for house escape matrix
	/// ordered members : 1) infectious, 2) susceptibles
	void init_house_for_evaluation(Household &house) override final
	{
		auto &order = cache_size_vector(house, HHsizevecKey::order);
		order.clear();
		order.reserve(house.size());
		order.insert(order.end(), house.infected_hhindexes().begin(), house.infected_hhindexes().end());
		order.insert(order.end(), house.susceptible_hhindexes().begin(), house.susceptible_hhindexes().end());
	}

	/// Risk of infection from the community.
	double community_risk() const { return par(ParamKey::alpha); }
	/// Risk of infection from another household member.
	double person_to_person_risk() const { return par(ParamKey::beta); }
	/// Store person-to-person escape probabilities inside a cache.
	void fill_house_escape_matrix(Household &house);

	/// Compute the log-likelihood for a single household.
	double compute_house_obs_log_lik(Household &house) override final;

	/// Simulate the infection within a household.
	void simulate_single_house(Household &house) override final;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Model : compute single house loglike  /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Model::compute_house_obs_log_lik(Household &house)
{
	// fill per-individual escape rate matrices
	fill_house_escape_matrix(house);

	const auto &order = cache_size_vector(house, HHsizevecKey::order);
	const auto &commu_mtx = cache_double_vector(house, HHdoublevecKey::commu_escape_rate);
	const auto &p2p_mtx = cache_double_matrix(house, HHdoublemtxKey::p2p_escape_rate);
	const unsigned n_infec = house.n_infected_status();
	const double commu_rate = community_risk();
	const double p2p_rate = person_to_person_risk();

	// identify distinct types among seropositives
	using type_key = std::tuple<double, double, double>;
	std::map<type_key, unsigned> type_map;
	std::vector<unsigned> type_counts;
	
	std::vector<double> commu_esc_rate_type;
	std::vector<double> infectivity_type;
	std::vector<double> susceptibility_type;

	for (unsigned i = 0; i < n_infec; ++i)
	{
		const Individual &member_i = house.member(order[i]);
		double com_i = member_i.group_rate(GroupRateKey::community);
		double sus_i = member_i.group_rate(GroupRateKey::susceptibility);
		double inf_i = member_i.group_rate(GroupRateKey::infectivity);

		type_key key = std::make_tuple(com_i, sus_i, inf_i);
		auto it = type_map.find(key);
		if (it == type_map.end())
		{
			unsigned tid = type_map.size();
			type_map[key] = tid;
			type_counts.push_back(1);
			commu_esc_rate_type.push_back(commu_rate * com_i);
			infectivity_type.push_back(inf_i);
			susceptibility_type.push_back(sus_i);
		}
		else
		{
			type_counts[it->second]++;
		}
	}

	// build K x K p2p escape rate matrix by type
	unsigned K = type_counts.size();
	std::vector<std::vector<double>> p2p_esc_rate_type(K, std::vector<double>(K));
	for (unsigned sigma = 0; sigma < K; ++sigma)
		for (unsigned tau = 0; tau < K; ++tau)
			p2p_esc_rate_type[sigma][tau] = p2p_rate
				* infectivity_type[sigma] * susceptibility_type[tau];

	return bino::chain_binom_grouped_log_lik(
		n_infec, type_counts,
		commu_esc_rate_type, p2p_esc_rate_type,
		commu_mtx, p2p_mtx);
}

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
		for (size_t k = 0; k < i; ++k)
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
