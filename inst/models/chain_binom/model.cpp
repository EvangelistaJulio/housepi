// model.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * CHAIN BINOMIAL * *** *  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//-----------------------------------------------------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

/* Community risk of infection (-logarithm of community escape probability). */
double Model::community_risk_of_infection() const
{
	return param("alpha");
}

/* Hazard of infection (-logarithm of person-to-person escape probability). */
double Model::hazard_of_infection() const
{
	return param("beta");
}

//------------------------------------------- /\ ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//-----------------------------------------------------------------------------------------------//

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Initialize the Lookup Table. */
void Model::init_model_for_inference()
{
	// Unknown status are not modeled
	if (report_data("n.miss.stat") > 0)
	{
		tools::error("A missing status was found : incompatible with the model.");
	}

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
double Model::compute_house_loglike(const Household &house) const
{
	int hh_size = house.size();
	int n_positive = house.num_infectious();

	double log_commu_esc_proba = -community_risk_of_infection();
	double log_commu_infec_proba = log(1. - exp(log_commu_esc_proba));
	double log_p2p_esc_proba = -hazard_of_infection();

	double log_likelihood = 0.;

	//================================================//
	//======== Contribution of Sero-Negatives ========//
	//================================================//

	log_likelihood += (hh_size - n_positive) * (log_commu_esc_proba + n_positive * log_p2p_esc_proba);

	//================================================//
	//======== Contribution of Sero-Positives ========//
	//================================================//

	if (n_positive > 0)
	{
		// sum over all possible generations
		double sero_positive_term = 0., tmp = 0.;
		for (unsigned c = 0; c < bax::Num_Gener_LookupTable[n_positive]; c++)
		{
			// reference to the LUT
			bax::gentype &gen_counts = bax::Generations_LookupTable[n_positive][c];

			// GENERATION (0) : infected by community
			tmp = gen_counts[0] * log_commu_infec_proba;
			// GENERATION (0) : escaped by community
			tmp += (n_positive - gen_counts[0]) * log_commu_esc_proba;

			// GENERATION (i) : within house, i > 0
			for (size_t i = 1; i < gen_counts.size(); i++)
			{
				// GENERATION (i) : infected by previous generation i - 1
				tmp += gen_counts[i] * log(1. - exp(gen_counts[i - 1] * log_p2p_esc_proba));

				// GENERATION (i) : escaped from other generations g, i. e. 0 <= g < i - 1
				for (size_t j = 0; j < i - 1; j++)
				{
					tmp += gen_counts[i] * gen_counts[j] * log_p2p_esc_proba;
				}
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
	double commu_esc_proba = exp(-community_risk_of_infection());
	double log_p2p_esc_proba = -hazard_of_infection();

	//==================================================//
	//===== GENERATION 0 : infection by community ======//
	//==================================================//

	for (int i = 0; i < house.size(); i++)
	{
		if (uniform() < commu_esc_proba)
		{
			house.change_status(i, SUSCEPTIBLE);
		}
		else
		{
			house.change_status(i, INFECTIOUS);
		}
	}

	//==============================================================//
	//===== GENERATION >=1 : infection by previous generation ======//
	//==============================================================//

	int num_old_gen = 0;
	int num_new_gen = house.num_infectious();

	while ((house.num_infectious() > 0) && (num_new_gen > 0))
	{
		num_old_gen = num_new_gen;
		num_new_gen = 0;

		// note : here househ[h].susceptible_set is not allowed because of change_status
		for (int i = 0; i < house.size(); i++)
		{
			if (house.status(i) == SUSCEPTIBLE)
			{
				if (uniform() >= exp(num_old_gen * log_p2p_esc_proba))
				{
					house.change_status(i, INFECTIOUS);
					num_new_gen++;
				}
			}
		}
	}
}