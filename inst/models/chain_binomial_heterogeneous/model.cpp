/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iterator>

#include "housepi_model_setup.h"
#include "binomial.h"
#include "model.h"

EMPTY_HOUSEPI_KEYS(AugmKey)
EMPTY_HOUSEPI_KEYS(DataKey)
EMPTY_HOUSEPI_KEYS(LongKey)
REGISTER_HOUSEPI_KEYS(GroupRateKey, community, susceptibility, infectivity)
EMPTY_HOUSEPI_KEYS(HHsizeKey)
REGISTER_HOUSEPI_KEYS(HHsizevecKey, order)
REGISTER_HOUSEPI_KEYS(HHdoublevecKey, commu_escape_rate)
REGISTER_HOUSEPI_KEYS(HHdoublemtxKey, p2p_escape_rate)
REGISTER_HOUSEPI_KEYS(ParamKey, alpha, beta)
HOUSEPI_MODEL_CONSTRUCTOR()

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : initialization  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::init_model_for_evaluation()
{
	// unknown status are not modeled
	assert_no_missing_status();

	// initialize look-up-table to largest possible value
	bax::init_escape_patterns_lookuptable(get_max_hh_cases());
}

void Model::init_house_for_evaluation(Household &house)
{
	auto &order = house.cache_size_vector(HHsizevecKey::order);
	order.clear();
	order.reserve(house.size());
	order.insert(order.end(), house.infected_hhindexes().begin(), house.infected_hhindexes().end());
	order.insert(order.end(), house.susceptible_hhindexes().begin(), house.susceptible_hhindexes().end());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Model : compute single house loglike  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Model::community_risk() const { return par(ParamKey::alpha); }

double Model::person_to_person_risk() const { return par(ParamKey::beta); }

double Model::compute_house_obs_log_lik(Household &house)
{
	// fill matrix with escape probabilities
	fill_house_escape_matrix(house);

	// compute corresponding loglikelihood times house multiplicity
	return house.multiplicity() * bax::chain_binomial_heter_log_lik(
									  house.n_infected_status(),
									  house.cache_double_vector(HHdoublevecKey::commu_escape_rate),
									  house.cache_double_matrix(HHdoublemtxKey::p2p_escape_rate));
}

void Model::fill_house_escape_matrix(Household &house) const
{
	const auto &order = house.cache_size_vector(HHsizevecKey::order);
	const double commu_rate = community_risk();
	const double p2p_rate = person_to_person_risk();

	// assign values to escape matrix
	auto &commu_mtx = house.cache_double_vector(HHdoublevecKey::commu_escape_rate);
	auto &p2p_mtx = house.cache_double_matrix(HHdoublemtxKey::p2p_escape_rate);

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
			// [i] : infector -> [r] : infectee
			p2p_mtx[i][k] = inf_i * p2p_rate * member_k.group_rate(GroupRateKey::susceptibility);
			// [r] : infector -> [i] : infectee
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
		if (uniform() < exp(-commu_rate * member.group_rate(GroupRateKey::community)))
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

		const auto sus_hhindexes = house.susceptible_hhindexes();
		for (size_t i : sus_hhindexes)
		{
			Individual &member_i = house.member(i);
			double sus_i = member_i.group_rate(GroupRateKey::susceptibility);

			double r = 0.;
			for (Individual *member_k : old_gen_set)
				r += member_k->group_rate(GroupRateKey::infectivity);

			if (uniform() >= exp(-r * p2p_rate * sus_i))
			{
				house.set_infected_status(member_i);
				new_gen_set.push_back(&member_i);
			}
		}
	}
}
