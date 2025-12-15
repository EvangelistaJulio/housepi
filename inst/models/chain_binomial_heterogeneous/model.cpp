/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 * 
 * Author: Francesco Camaglia
 * Copyright © 2025 Institut Pasteur (Paris), CNRS and INSERM.
 * See the COPYRIGHT file for details.
 * 
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the COPYING file for details.
 *
 * model.cpp
 * Model template of the 'housepi' package.
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"
#include "binomial.h"

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
	bax::init_escape_patterns_lookuptable(get_report_hh_size().max_cases.first);

	// allocate memory for the escape matrix
	cache_house_size_t_vector("order");
	cache_house_double_vector("commu-escape-proba");
	cache_house_double_matrix("p2p-escape-proba");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Model : compute single house loglike  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Compute the log-likelihood for a single household.
double Model::compute_household_obs_log_lik(Household &house)
{
	// fill matrix with escape probabilities
	fill_house_escape_matrix(house);

	// compute corresponding loglikelihood times house multiplicity
	return house.multiplicity() * bax::chain_binomial_heter_log_lik(
									  house.n_infec_status(),
									  house.cache_double_vector("commu-escape-proba"),
									  house.cache_double_matrix("p2p-escape-proba"));
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

	// ordered members group rates : 1) infectious, 2) susceptibles
	auto &order = house.cache_size_t_vector("order");
	// note 1: this works only if size = n_infec_status + n_healthy_status
	// note 2: this is not changing throughout inference w\o status change
	order = house.infected_hhindexes() + house.healthy_hhindexes();

	// assign values to escpae matrix
	auto &commu_mtx = house.cache_double_vector("commu-escape-proba");
	auto &p2p_mtx = house.cache_double_matrix("p2p-escape-proba");

	// fill matrix according to order
	for (size_t i = 0; i < house.size(); ++i)
	{
		const Individual &member_i = house.get_member(order[i]);
		double inf_i = member_i.group_rate("infectivity");
		double sus_i = member_i.group_rate("susceptibility");

		// community
		commu_mtx[i] = std::exp(-commu_rate * member_i.group_rate("community"));
		// infector always "escapes" himself (check if actually relevant)
		p2p_mtx[i][i] = 1.0;

		// person-to-person
		for (size_t k = 0; k < i; ++k)
		{
			const Individual &member_k = house.get_member(order[k]);
			// [i] : infector -> [r] : infectee
			p2p_mtx[i][k] = std::exp(-inf_i * p2p_rate * member_k.group_rate("susceptibility"));
			// [r] : infector -> [i] : infectee
			p2p_mtx[k][i] = std::exp(-member_k.group_rate("infectivity") * p2p_rate * sus_i);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////// Model : transmission rate  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Risk of infection from the community.
double Model::community_risk() const
{
	return par("alpha");
}

/// Risk of infection from another household member.
double Model::person_to_person_risk() const
{
	return par("beta");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/// Simulate the infection within a household.
void Model::simulate_single_household(Household &house)
{
	const double commu_rate = community_risk();
	const double p2p_rate = person_to_person_risk();

	//////////////////////////////////////////////////////////////////
	/////////////  Generation 0 : infection by community  ////////////
	//////////////////////////////////////////////////////////////////

	std::vector<Individual *> new_gen_set = {};

	for (Individual &member : house)
		if (uniform() < exp(-commu_rate * member.group_rate("community")))
			house.set_healthy_status(member);
		else
		{
			house.set_infected_status(member);
			new_gen_set.push_back(&member);
		}

	//////////////////////////////////////////////////////////////////
	//////  Generation >=1 : infection by previous generation  ///////
	//////////////////////////////////////////////////////////////////

	std::vector<Individual *> old_gen_set = {};

	while ((house.n_healthy_status() > 0) && (!new_gen_set.empty()))
	{
		std::swap(old_gen_set, new_gen_set);
		new_gen_set.clear();

		for (Individual &member_i : house.healthy())
		{
			double sus_i = member_i.group_rate("susceptibility");

			double r = 0.;
			for (Individual *member_k : old_gen_set)
				r += member_k->group_rate("infectivity");

			if (uniform() >= exp(-r * p2p_rate * sus_i))
			{
				house.set_infected_status(member_i);
				new_gen_set.push_back(&member_i);
			}
		}
	}
}
