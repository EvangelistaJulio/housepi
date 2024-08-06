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
double Model::community_risk_of_infection()
{
	return get_value("alpha");
}

/* Hazard of infection (-logarithm of person-to-person escape probability). */
double Model::hazard_of_infection()
{
	return get_value("beta");
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
	if (report_data["n.miss.stat"] > 0)
	{
		screen::raise_error("A missing status was found : incompatible with the model.");
	}

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_generations_lookuptable(report_data["max.n.HH.cases"]);
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
double Model::compute_single_house_loglike(int h)
{
	double log_likelihood = 0;

	int hh_size = househ[h].house_size;
	int n_positive = househ[h].num_infectious;

	double log_commu_esc_proba = -community_risk_of_infection();
	double log_commu_infec_proba = log(1. - exp(log_commu_esc_proba));
	double log_p2p_esc_proba = -hazard_of_infection();

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
		double  sero_positive_term = 0., tmp = 0.;
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
void Model::simulate_single_house(int h)
{
	double commu_esc_proba = exp(-community_risk_of_infection());
	double log_p2p_esc_proba = -hazard_of_infection();

	//================================================//
	//==== GENERATION 0 : infection by community =====//
	//================================================//

	for (int i = 0; i < househ[h].house_size; i++)
	{
		if (uniform() <= commu_esc_proba)
		{
			househ[h].change_status(i, SUSCEPTIBLE);
		}		
		else
		{
			househ[h].change_status(i, INFECTIOUS);
		}
	}

	//============================================================//
	//==== GENERATION >=1 : infection by previous generation =====//
	//============================================================//

	int num_gen_infec = househ[h].num_infectious;
	while ((househ[h].num_susceptible > 0) && (num_gen_infec > 0))
	{
		// note : here househ[h].susceptible_set is not allowed because of change_status
		for (int i = 0; i < househ[h].house_size; i++)
		{
			if ((househ[h].is_susceptible[i]) && (uniform() > exp(num_gen_infec * log_p2p_esc_proba)))
			{
				househ[h].change_status(i, INFECTIOUS);
			}
		}
		num_gen_infec = househ[h].num_infectious - num_gen_infec;
	}
}