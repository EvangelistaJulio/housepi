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

#include "housepi.h"
#include "binomial.h"
#include "mcmc.h"

HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE()
HOUSEPI_DECLARE_EMPTY_KEYS(AugmKey)
HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
HOUSEPI_DECLARE_REGISTER_KEYS(GroupRateKey, community, susceptibility, infectivity)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizeKey)
HOUSEPI_DECLARE_REGISTER_KEYS(HHsizevecKey, order)
HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublevecKey, commu_escape_rate)
HOUSEPI_DECLARE_REGISTER_KEYS(HHdoublemtxKey, p2p_escape_rate)
HOUSEPI_DECLARE_REGISTER_KEYS(ParamKey, alpha, beta)

class Model final : public MCMC
{
public:
	HOUSEPI_MODEL_CONSTRUCTOR()

	/// Initialize the LookUp Tables to largest possible values.
	/// NOTE: missing status are not modeled.
	void init_model_for_evaluation() override final
	{
		assert_no_missing_status();
		bax::init_escape_patterns_lookuptable(get_max_hh_cases());
	}

	/// Initialize internal order for house escape matrix
	/// ordered members : 1) infectious, 2) susceptibles
	/// NOTE: this works only if size = n_infected_status + n_susceptible_status
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
	/// Important: it must be re-ordered to have seropositives first.
	void fill_house_escape_matrix(Household &house);

	/// Compute the log-likelihood for a single household.
	double compute_house_obs_log_lik(Household &house) override final;

	/// Simulate the infection within a household.
	void simulate_single_house(Household &house) override final;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Model : compute single house loglike  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Model::compute_house_obs_log_lik(Household &house)
{
	// fill matrix with escape probabilities
	fill_house_escape_matrix(house);

	// compute corresponding loglikelihood times house multiplicity
	return bax::chain_binomial_heter_log_lik(
		house.n_infected_status(),
		cache_double_vector(house, HHdoublevecKey::commu_escape_rate),
		cache_double_matrix(house, HHdoublemtxKey::p2p_escape_rate));
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

HOUSEPI_CLOSE_NAMED_MODULE_NAMESPACE()
HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE()
