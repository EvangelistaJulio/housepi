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

#include "housepi_model_setup.h"
#include "binomial.h"
#include "model.h"

HOUSEPI_DECLARE_EMPTY_KEYS(AugmKey)
HOUSEPI_DECLARE_EMPTY_KEYS(DataKey)
HOUSEPI_DECLARE_EMPTY_KEYS(LongKey)
HOUSEPI_DECLARE_EMPTY_KEYS(GroupRateKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizeKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHsizevecKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublevecKey)
HOUSEPI_DECLARE_EMPTY_KEYS(HHdoublemtxKey)
HOUSEPI_DECLARE_REGISTER_KEYS(ParamKey, alpha, beta)
HOUSEPI_DEFINE_MODEL_CONSTRUCTOR()
HOUSEPI_DECLARE_RCPP_MODULE(HOUSEPI_MODULE_NAME)

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////  Model : initialization  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::init_model_for_evaluation()
{
	// missing status are not modeled
	assert_no_missing_status();

	// initialize look-up-table to largest possible value
	bax::init_generations_lookuptable(get_max_hh_cases());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Model::community_risk() const { return par(ParamKey::alpha); }

double Model::person_to_person_risk() const { return par(ParamKey::beta); }

double Model::compute_house_obs_log_lik(Household &house)
{
	double log_likelihood = 0.;

	double log_commu_esc_proba = -community_risk();
	double log_commu_infec_proba = log(1. - exp(log_commu_esc_proba));
	double log_p2p_esc_proba = -person_to_person_risk();

	auto n_positive = house.n_infected_status();

	///////////////////////////////////////////
	/////////  Sero-Negatives Members  ////////
	///////////////////////////////////////////

	log_likelihood += (house.size() - n_positive) * (log_commu_esc_proba + n_positive * log_p2p_esc_proba);

	///////////////////////////////////////////
	/////////  Sero-Positives Members  ////////
	///////////////////////////////////////////

	if (n_positive > 0)
	{
		// sum over all possible generations
		double sero_positive_term = 0., tmp = 0.;
		for (unsigned long c = 0; c < bax::get_num_generations(n_positive); ++c)
		{
			// reference to the LUT
			const auto &gen_counts = bax::get_generation_vector(n_positive, c);

			// generation (0) : infected by community
			tmp = gen_counts[0] * log_commu_infec_proba;
			// generation (0) : escaped by community
			tmp += (n_positive - gen_counts[0]) * log_commu_esc_proba;

			// generation (i) : within house, i > 0
			for (size_t i = 1; i < gen_counts.size(); ++i)
			{
				// generation (i) : infected by previous generation i - 1
				tmp += gen_counts[i] * std::log(1. - std::exp(gen_counts[i - 1] * log_p2p_esc_proba));

				// generation (i) : escaped from other generations g, i. e. 0 <= g < i - 1
				for (size_t j = 0; j < i - 1; ++j)
					tmp += gen_counts[i] * gen_counts[j] * log_p2p_esc_proba;
			}

			// times n. of equivalent configurations returning the same kind of generation
			sero_positive_term += bax::get_generation_multipl(n_positive, c) * std::exp(tmp);
		}

		log_likelihood += std::log(sero_positive_term);
	}

	return log_likelihood;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : simulate single house  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::simulate_single_house(Household &house)
{
	double commu_esc_proba = exp(-community_risk());
	double log_p2p_esc_proba = -person_to_person_risk();

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

	std::vector<size_t> sus_hhindexes;
	sus_hhindexes.reserve(house.n_susceptible_status());
	while ((house.n_infected_status() > 0) && (n_new_gen > 0))
	{
		n_old_gen = n_new_gen;
		n_new_gen = 0;

		sus_hhindexes = house.susceptible_hhindexes();
		for (size_t i : sus_hhindexes)
			if (uniform() >= exp(n_old_gen * log_p2p_esc_proba))
			{
				house.set_infected_status(house.member(i));
				++n_new_gen;
			}
	}
}
