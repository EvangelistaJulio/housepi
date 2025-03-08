/*
 * model.cpp
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  * *** * REED-FROST SINGLE INTRODUCTION * *** *  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//-----------------------------------------------------------------------------------------------//
//---------------------------------------- CHANGE-ME --------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

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

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Model : common initialization  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the Lookup Table. */
void Model::init_model_for_inference()
{
	// Unknown status are not modeled
	if (get_report_indiv("missing") > 0)
		tools::error("A missing status was found : incompatible with the model.");

	// households without case don't have a contribution
	if (get_report_house("no cases") == 0)
		tools::warning(
			"Hoseholds without cases are not modeled and will be ignored. Discard them to improve computation time.");

	// households of size 1 don't have a contribution
	if (get_report_hh_size("min size").first == 1)
		tools::warning(
			"Hoseholds with one single member are not modeled and will be ignored. Discard them to improve computation time.");

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_generations_lookuptable(get_report_hh_size("max cases").first);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Model : Methods (likelihood)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(Household &house) const
{
	double log_likelihood = 0;
	double log_p2p_esc_proba = -person_to_person_risk();
	double log_p2p_infec_proba = log(1. - exp(log_p2p_esc_proba));

	indiv n_positive = house.n_infec_status();

	//================================================//
	//======== Contribution of Sero-Negatives ========//
	//================================================//

	log_likelihood += (house.size() - n_positive) * n_positive * log_p2p_esc_proba;

	//================================================//
	//======== Contribution of Sero-Positives ========//
	//================================================//

	if (n_positive > 1) // single introduction
	{
		// sum over all possible generations
		double sero_positive_term = 0., tmp = 0.;
		for (unsigned c = 0; c < bax::Num_Gener_LookupTable[n_positive - 1]; ++c)
		{
			// reference to the LUT
			bax::gentype &gen_counts = bax::Generations_LookupTable[n_positive - 1][c];

			// GENERATION (0) : infected by case zero
			tmp = gen_counts[0] * log_p2p_infec_proba;
			// GENERATION (0) : escaped by case zero
			tmp += (n_positive - 1 - gen_counts[0]) * log_p2p_esc_proba;

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
			sero_positive_term += bax::Multipl_Gener_LookupTable[n_positive - 1][c] * exp(tmp);
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
///////////////////////////////  Model : Methods (simulation)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Simulate the infection within a household. */
void Model::simulate_single_house(Household &house)
{
	double log_p2p_esc_proba = -person_to_person_risk();

	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	// set all cases to susceptible
	for (indiv i = 0; i < house.size(); i++)
		house.change_status(i, HEALTHY_STATUS);

	// pick a random introduction
	indiv first_case = uniform() * house.size();
	house.change_status(first_case, INFECTED_STATUS);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	ushort n_old_gen = 0;
	ushort n_new_gen = 1;

	while ((house.n_healthy_status() > 0) && (n_new_gen > 0))
	{
		n_old_gen = n_new_gen;
		n_new_gen = 0;

		// note : here house.healthy_status_list is not allowed because of change_status
		for (indiv i = 0; i < house.size(); ++i)
			if (house.status(i) == HEALTHY_STATUS)
				if (uniform() >= exp(n_old_gen * log_p2p_esc_proba))
				{
					house.change_status(i, INFECTED_STATUS);
					++n_new_gen;
				}
	}
}
