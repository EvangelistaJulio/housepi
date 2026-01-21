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

REGISTER_HOUSEPI_KEYS(AugmKey, )
REGISTER_HOUSEPI_KEYS(DataKey, )
REGISTER_HOUSEPI_KEYS(GroupRateKey, community, susceptibility, infectivity)
REGISTER_HOUSEPI_KEYS(HHsize_tKey, order)
REGISTER_HOUSEPI_KEYS(HHdoubleKey, commu_escape_proba)
REGISTER_HOUSEPI_KEYS(HHmatrixKey, p2p_escape_proba)
REGISTER_HOUSEPI_KEYS(ParamKey, alpha, beta)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : initialization  /////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Initialize the Lookup Table.
void Model::init_model_for_inference()
{
	// unknown status are not modeled
	assert_no_missing_status();

	// initialize look-up-table to largest possible value
	bax::init_escape_patterns_lookuptable(get_max_hh_cases());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Model : compute single house loglike  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute the log-likelihood for a single household.
double Model::compute_house_obs_log_lik(Household &house)
{
	// fill matrix with escape probabilities
	fill_house_escape_matrix(house);

	// compute corresponding loglikelihood times house multiplicity
	return house.multiplicity() * bax::chain_binomial_heter_log_lik(
									  house.n_infected_status(),
									  house.cache_double_vector(HHdoubleKey::commu_escape_proba),
									  house.cache_double_matrix(HHmatrixKey::p2p_escape_proba));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Model : fill house escape matrix  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// The person-to-person escape probabilities stored inside a matrix.
/// Important: it must be re-ordered to have seropositives first.
void Model::fill_house_escape_matrix(Household &house) const
{
	const double commu_rate = community_risk();
	const double p2p_rate = person_to_person_risk();

	/// ordered members group rates : 1) infectious, 2) susceptibles
	/// NOTE: this works only if size = n_infected_status + n_susceptible_status
	/// UPDATE: this is not changing throughout inference w\o status change!
	auto &order = house.cache_size_t_vector(HHsize_tKey::order);
	order.reserve(house.size());
	order.insert(order.end(), house.infected_hhindexes().begin(), house.infected_hhindexes().end());
	order.insert(order.end(), house.healthy_hhindexes().begin(), house.healthy_hhindexes().end());

	// assign values to escape matrix
	auto &commu_mtx = house.cache_double_vector(HHdoubleKey::commu_escape_proba);
	auto &p2p_mtx = house.cache_double_matrix(HHmatrixKey::p2p_escape_proba);

	// fill matrix according to order
	for (size_t i = 0; i < house.size(); ++i)
	{
		const auto &member_i = house.get_member(order[i]);
		double inf_i = member_i.group_rate(GroupRateKey::infectivity);
		double sus_i = member_i.group_rate(GroupRateKey::susceptibility);

		// community
		commu_mtx[i] = std::exp(-commu_rate * member_i.group_rate(GroupRateKey::community));
		// infector always "escapes" himself (check if actually relevant)
		p2p_mtx[i][i] = 1.0;

		// person-to-person
		for (size_t k = 0; k < i; ++k)
		{
			const auto &member_k = house.get_member(order[k]);
			// [i] : infector -> [r] : infectee
			p2p_mtx[i][k] = std::exp(-inf_i * p2p_rate * member_k.group_rate(GroupRateKey::susceptibility));
			// [r] : infector -> [i] : infectee
			p2p_mtx[k][i] = std::exp(-member_k.group_rate(GroupRateKey::infectivity) * p2p_rate * sus_i);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Model : transmission rate  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Risk of infection from the community.
double Model::community_risk() const { return par(ParamKey::alpha); }

/// Risk of infection from another household member.
double Model::person_to_person_risk() const { return par(ParamKey::beta); }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Simulate the infection within a household.
void Model::simulate_single_house(Household &house)
{
	const double commu_rate = community_risk();
	const double p2p_rate = person_to_person_risk();

	//////////////////////////////////////////////////////////////////
	/////////////  Generation 0 : infection by community  ////////////
	//////////////////////////////////////////////////////////////////

	std::vector<Individual *> new_gen_set = {};

	for (auto &member : house)
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

		for (auto &member_i : house.healthy())
		{
			double sus_i = member_i.group_rate(GroupRateKey::susceptibility);

			double r = 0.;
			for (auto *member_k : old_gen_set)
				r += member_k->group_rate(GroupRateKey::infectivity);

			if (uniform() >= exp(-r * p2p_rate * sus_i))
			{
				house.set_infected_status(member_i);
				new_gen_set.push_back(&member_i);
			}
		}
	}
}
