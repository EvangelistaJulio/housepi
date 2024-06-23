// model.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  * *** * CHAIN BINOMIAL * *** *  ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//------------------------------------------------------------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- \/ -------------------------------------------------------//

/* The logarithm of the community escape probability. */
double Model::log_community_escape_probability()
{
	return -get_value("alpha");
}

/* The logarithm of the person-to-person escape probability. */
double Model::log_person_to_person_escape_probability()
{
	return -get_value("beta");
}

//------------------------------------------------- /\ -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/* Initialize the Lookup Table. */
void Model::initialize_model()
{
	// Unknown status are not modeled
	if (report_data["n.unk.stat"] > 0)
	{
		raise_error_message("A missing status was found : incompatible with the model.");
	}

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_generations_lookuptable(report_data["max.n.HH.cases"]);
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(int h)
{
	double log_likelihood = 0;

	int hh_size = househ[h].house_size;
	int n_positive = househ[h].infectious_list.size();

	double log_commu_esc_proba = log_community_escape_probability();
	double log_commu_infec_proba = log(1. - exp(log_commu_esc_proba));
	double log_p2p_esc_proba = log_person_to_person_escape_probability();

	//==============================================//
	//========Contribution of Sero-Negatives========//
	//==============================================//

	log_likelihood += (hh_size - n_positive) * (log_commu_esc_proba + n_positive * log_p2p_esc_proba);

	//==============================================//
	//========Contribution of Sero-Positives========//
	//==============================================//

	if (n_positive > 0)
	{
		int counter = 0;
		double tmp = 0., sero_positive_term = 0.;

		// sum over all possible generations
		for (bax::gentype gen_counts : bax::Generations_LookupTable[n_positive])
		{
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
			sero_positive_term += bax::Num_Generations_LookupTable[n_positive][counter] * exp(tmp);
			counter++;
		}

		log_likelihood += log(sero_positive_term);
	}

	return log_likelihood;
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////  Model : Methods (simulation)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Simulate the infection within a household. */
void Model::simulate_single_house(int h)
{
	int num_new_infectious = 0;
	int num_old_infectious = 0;
	std::vector<int> store_old_list = {};

	double commu_esc_proba = exp(log_community_escape_probability());
	double log_p2p_esc_proba = log_person_to_person_escape_probability();

	//====================//
	//=====INITIALIZE=====//
	//====================//

	// set all cases to healthy (change infectious and unknown to susceptibles)
	store_old_list = househ[h].infectious_list;
	for (int i : store_old_list)
	{
		househ[h].change_status(i, SUSCEPTIBLE);
	}

	//==============================================//
	//====GENERATION 0 : infection by community=====//
	//==============================================//

	store_old_list = househ[h].susceptibles_list;
	for (int i : store_old_list)
	{
		if (uniform() > commu_esc_proba)
		{
			househ[h].change_status(i, INFECTIOUS);
			num_new_infectious++;
		}
	}
	num_new_infectious = 0;

	//==========================================================//
	//====GENERATION >=1 : infection by previous generation=====//
	//==========================================================//

	while ((househ[h].susceptibles_list.size() > 0) && (num_new_infectious > 0))
	{
		num_old_infectious = num_new_infectious;
		num_new_infectious = 0;

		store_old_list = househ[h].susceptibles_list;
		for (int i : store_old_list)
		{
			if (uniform() > exp(num_old_infectious * log_p2p_esc_proba))
			{
				househ[h].change_status(i, INFECTIOUS);
				num_new_infectious++;
			}
		}
	}
}