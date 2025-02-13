/*
 * model.cpp
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * CHAIN BINOMIAL * *** *  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//-----------------------------------------------------------------------------------------------//
//---------------------------------------- CHANGE-ME --------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

/* [CHANGE-ME] Risk of infection from the community. */
double Model::community_risk() const
{
	return param("alpha");
}

/* [CHANGE-ME] Risk of infection from another household member. */
double Model::person_to_person_risk() const
{
	return param("beta");
}

//------------------------------------------- /\ ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//---------------------------------------- CHANGE-ME --------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Initialize the Lookup Table. */
void Model::init_model_for_inference()
{
	// Unknown status are not modeled
	if (report_data("n.miss.stat") > 0)
		tools::error("A missing status was found : incompatible with the model.");

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_generations_lookuptable(report_data("max.n.HH.cases"));
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : Methods (likelihood)  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(Household &house) const
{
	double log_likelihood = 0.;

	double log_commu_esc_proba = -community_risk();
	double log_commu_infec_proba = log(1. - exp(log_commu_esc_proba));
	double log_p2p_esc_proba = -person_to_person_risk();

	ushort n_positive = house.n_infectious();

	//================================================//
	//======== Contribution of Sero-Negatives ========//
	//================================================//

	log_likelihood += (house.size() - n_positive) * (log_commu_esc_proba + n_positive * log_p2p_esc_proba);

	//================================================//
	//======== Contribution of Sero-Positives ========//
	//================================================//

	if (n_positive > 0)
	{
		// sum over all possible generations
		double sero_positive_term = 0., tmp = 0.;
		for (unsigned c = 0; c < bax::Num_Gener_LookupTable[n_positive]; ++c)
		{
			// reference to the LUT
			bax::gentype &gen_counts = bax::Generations_LookupTable[n_positive][c];

			// GENERATION (0) : infected by community
			tmp = gen_counts[0] * log_commu_infec_proba;
			// GENERATION (0) : escaped by community
			tmp += (n_positive - gen_counts[0]) * log_commu_esc_proba;

			// GENERATION (i) : within house, i > 0
			for (size_t i = 1; i < gen_counts.size(); ++i)
			{
				// GENERATION (i) : infected by previous generation i - 1
				tmp += gen_counts[i] * log(1. - exp(gen_counts[i - 1] * log_p2p_esc_proba));

				// GENERATION (i) : escaped from other generations g, i. e. 0 <= g < i - 1
				for (size_t j = 0; j < i - 1; ++j)
					tmp += gen_counts[i] * gen_counts[j] * log_p2p_esc_proba;
			}

			// times n. of equivalent configurations returning the same kind of generation
			sero_positive_term += bax::Multipl_Gener_LookupTable[n_positive][c] * exp(tmp);
		}

		log_likelihood += log(sero_positive_term);
	}

	return log_likelihood;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : Methods (simulation)  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Simulate the infection within a household. */
void Model::simulate_single_house(Household &house)
{
	double commu_esc_proba = exp(-community_risk());
	double log_p2p_esc_proba = -person_to_person_risk();

	//==================================================//
	//===== GENERATION 0 : infection by community ======//
	//==================================================//

	for (indiv i = 0; i < house.size(); ++i)
		if (uniform() < commu_esc_proba)
			house.change_status(i, HEALTHY_STATUS);
		else
			house.change_status(i, INFECTED_STATUS);

	//==============================================================//
	//===== GENERATION >=1 : infection by previous generation ======//
	//==============================================================//

	ushort n_old_gen = 0;
	ushort n_new_gen = house.n_infectious();

	while ((house.n_infectious() > 0) && (n_new_gen > 0))
	{
		n_old_gen = n_new_gen;
		n_new_gen = 0;

		// note : loop over 'house.susceptible_set()' is not allowed because of 'change_status'
		for (indiv i = 0; i < house.size(); ++i)
		{
			if (house.status(i) == HEALTHY_STATUS)
			{
				if (uniform() >= exp(n_old_gen * log_p2p_esc_proba))
				{
					house.change_status(i, INFECTED_STATUS);
					++n_new_gen;
				}
			}
		}
	}
}
