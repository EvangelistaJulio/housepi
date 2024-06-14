// model.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  * *** * CHAIN BINOMIAL with SINGLE INTRODUCTION * *** *  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//------------------------------------------------------------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- \/ -------------------------------------------------------//

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

// auxiliary global variable
std::map<int, std::vector<gentype>> Generations_Lookup_Table = {};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////  WARNINGS  /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

bool Single_Introduction_Warning = false;
bool Single_Member_House_Warning = false;

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// Auxiliary Functions ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Standard definition for the factorial. */
unsigned factorial(unsigned n)
{
	unsigned result = 1;
	for (unsigned i = 2; i <= n; i++)
	{
		result *= i;
	}
	return result;
}

/* All possible ranked offspring generations for a given number of infected individuals n > 0.

	e.g.
	----
	n_infec = 3
	return [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 2]]
*/
std::vector<gentype> generation_span(int n)
{
	int root = n - 1;
	std::vector<gentype> output(1, gentype(n, 0));

	while (root > 0)
	{
		for (gentype g : generation_span(n - root))
		{
			gentype temp(root, 0);
			for (int i = 0; i < n - root; i++)
			{
				g[i] += 1;
			}
			temp.insert(temp.end(), g.begin(), g.end());
			output.push_back(temp);
		}
		root--;
	}

	return output;
}

/* Count-compression of all possible ranked offpsring generations. */
std::vector<gentype> generation_count_span(int n)
{
	std::vector<gentype> output = {};
	for (gentype g : generation_span(n))
	{
		// I take advantage of the order granted by the map
		std::map<int, int> freq_map = {};
		for (int val : g)
		{
			if (freq_map.find(val) != freq_map.end())
			{
				freq_map[val]++;
			}
			else
			{
				freq_map.emplace(val, 1);
			}
		}

		gentype counts = {};
		for (const auto &pair : freq_map)
		{
			counts.push_back(pair.second);
		}
		output.push_back(counts);
	}

	return output;
}

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

	// households without case don't have a contribution
	if (report_data["no.HH.case"] == 0)
	{
		raise_warning_message(
			"Hoseholds without cases are not modeled and will be ignored. Discard them to increase computation time.",
			Single_Introduction_Warning);
	}

	// households of size 1 don't have a contribution
	if (report_data["min.HH.size"] == 1)
	{
		raise_warning_message(
			"Hoseholds with one single member are not modeled and will be ignored. Discard them to increase computation time.", 
			Single_Member_House_Warning);
	}

	//=====Find maximum possible number of sero-positives=====//
	int max_num = 0, aux = 0;
	for (int h = 0; h < num_househ; h++)
	{
		aux = househ[h].infectious_list.size();
		max_num = (aux > max_num) ? aux : max_num;
	}

	//=====Initialize look-up-table to largest possible value=====//
	Generations_Lookup_Table.emplace(0, std::vector<gentype>());
	for (int n = 1; n <= max_num; n++)
	{
		Generations_Lookup_Table.emplace(n, generation_count_span(n));
	}
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

	double log_p2p_esc_proba = log_person_to_person_escape_probability();
	double log_p2p_infec_proba = log(1. - exp(log_p2p_esc_proba));

	//==============================================//
	//========Contribution of Sero-Negatives========//
	//==============================================//

	log_likelihood += (hh_size - n_positive) * n_positive * log_p2p_esc_proba;

	//==============================================//
	//========Contribution of Sero-Positives========//
	//==============================================//

	int multiplicity;
	double tmp, sero_positive_term;

	if (n_positive > 0)
	{
		sero_positive_term = 0.;

		// sum over all possible generations
		for (gentype gen_counts : Generations_Lookup_Table.at(n_positive - 1)) // remove single introduction
		{
			// number of equivalent configurations returning the same kind of generation
			multiplicity = factorial(n_positive) / factorial(gen_counts[0]); // euclidean division

			// GENERATION (0) : infected by community
			tmp = gen_counts[0] * log_p2p_infec_proba;
			// GENERATION (0) : escaped by community
			tmp += (n_positive - 1 - gen_counts[0]) * log_p2p_infec_proba;

			// GENERATION (i) : within house, i > 0
			for (size_t i = 1; i < gen_counts.size(); i++)
			{
				// GENERATION (i) : infected by previous generation i - 1
				tmp += gen_counts[i] * log(1. - exp(gen_counts[i - 1] * log_p2p_esc_proba));

				// GENERATION (i) : escaped from all others generation s.t. 0 <= generation < i - 1
				for (size_t j = 0; j < i - 1; j++)
				{
					tmp += gen_counts[i] * gen_counts[j] * log_p2p_esc_proba;
				}

				multiplicity /= factorial(gen_counts[i]);
			}

			sero_positive_term += multiplicity * exp(tmp);
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
	int hh_size = househ[h].house_size;

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
	//====GENERATION 0 : by single introduction=====//
	//==============================================//

	// pick a random introduction
	househ[h].change_status(static_cast<int>(uniform() * hh_size), INFECTIOUS);
	num_new_infectious = 1;

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
