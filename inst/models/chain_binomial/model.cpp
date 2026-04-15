/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
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
#define _HOUSEPI_MODULE chain_binomial
#endif

HOUSEPI_OPEN_MODULE_NAMESPACE()
class Model final : public ModelBase
{
public:
	HOUSEPI_DECLARE_EMPTY_KEYS(AugmKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(GroupRateKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHsizeKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHsizevecKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
	HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtxKey)
	HOUSEPI_DECLARE_REGISTER_KEYS(ParamKey, alpha, beta, gamma, size_ref)

	/// General advanced inscrition for the setup.
	void init_advanced_setup() override final
	{
		request_discrete_param(ParamKey::size_ref);
	}

	/// Initialize the LookUp Tables to largest possible values.
	/// NOTE: missing status are not modeled
	void init_model_for_evaluation() override final
	{
		assert_no_missing_status();
		bino::init_generations_lookuptable(get_max_hh_cases());
	}

	/// Risk of infection from the community.
	double community_risk() const { return par(ParamKey::alpha); }
	/// Risk of infection from another household member.
	double person_to_person_risk(size_t house_size) const
	{
		return par(ParamKey::beta) *
			   std::exp(par(ParamKey::gamma) * std::log(par(ParamKey::size_ref) / house_size));
	}

	/// Compute the log-likelihood within household.
	double compute_house_obs_log_lik(Household &house) override final
	{
		size_t house_size = house.size();
		return housepi::bino::chain_binom_homogeneous_log_lik(
			house.n_infected_status(), house_size,
			community_risk(), person_to_person_risk(house_size));
	}

	/// Simulate the infection within a household.
	void simulate_single_house(Household &house) override final;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::simulate_single_house(Household &house)
{
	double commu_esc_proba = std::exp(-community_risk());
	double p2p_rate = person_to_person_risk(house.size());

	//////////////////////////////////////////////////////////////////
	/////////////  Generation 0 : infection by community  ////////////
	//////////////////////////////////////////////////////////////////

	for (auto &member : house)
		if (uniform() < commu_esc_proba)
			house.set_susceptible_status(member);
		else
			house.set_infected_status(member);

	//////////////////////////////////////////////////////////////////
	//////  Generation >=1 : infection by previous generation  ///////
	//////////////////////////////////////////////////////////////////

	size_t n_old_gen = 0;
	size_t n_new_gen = house.n_infected_status();

	while (house.n_infected_status() > 0 && n_new_gen > 0)
	{
		n_old_gen = n_new_gen;
		n_new_gen = 0;

		for (Individual &member_i : house)
		{
			if (member_i.is_infected_status())
				continue;

			if (exponential() < n_old_gen * p2p_rate)
			{
				house.set_infected_status(member_i);
				++n_new_gen;
			}
		}
	}
}
HOUSEPI_CLOSE_MODULE_NAMESPACE()
