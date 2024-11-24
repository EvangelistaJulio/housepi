/*
 * model.cpp
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 *
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"
#include "binomial.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  Global Variables  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//----- Look-Up-Tables -----//
std::vector<std::vector<double>> LUT_commu_esc_proba;
std::vector<squarematrix> LUT_p2p_esc_proba;

//-----------------------------------------------------------------------------------------------//
//---------------------------------------- CHANGE-ME --------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

/* [CHANGE-ME] Assign the parameters controlling infectivity for each individual. */
void Model::assign_group_rate_to_individual(Household &house, int i)
{
	bool is_child = house.data(i, "age") <= param("child_age");

	//=====Susceptibility=====//
	if (is_child)
	{ // child
		house.assign_group_rate(i, "Susceptibility", "sus_child");
	}
	else
	{ // adult
		if (house.data(i, "sex") == 0)
		{ // adult women
			house.assign_group_rate(i, "Susceptibility", "sus_women");
		}
		else
		{ // adult men : Susceptibility [REFERENCE]
			house.assign_group_rate(i, "Susceptibility", REFERENCE);
		}
	}

	//=====Infectivity=====//
	if (!is_child)
	{ // adult : Susceptiblity [REFERENCE]
		house.assign_group_rate(i, "Infectivity", REFERENCE);
	}
	else
	{ // child
		house.assign_group_rate(i, "Infectivity", "inf_child");
	}
}

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

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Model : fill house escape matrix  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The person-to-person escape probabilities stored inside a matrix.
	Important: it must be re-ordered to have seropositives first. */
void Model::fill_house_escape_matrix(const Household &house) const
{
	double commu_transm_rate = community_risk();
	double p2p_transm_rate = person_to_person_risk();

	//----- ordered members group rates : 1) infectious, 2) susceptibles -----//
	std::vector<int> order;
	order.reserve(house.size());
	// 1) Append infectious members
	order.insert(order.end(), house.infectious_set().begin(), house.infectious_set().end());
	// 2) Append susceptible members
	order.insert(order.end(), house.susceptible_set().begin(), house.susceptible_set().end());

	//----- Assign values to LUT matrix -----//
	size_t H = house.index();
	for (int c = 0; c < house.size(); ++c)
	{
		double susc_commu_c = house.group_rate(order[c], "Susc-Community");
		double infect_c = house.group_rate(order[c], "Infectivity");
		double susc_house_c = house.group_rate(order[c], "Susc-House");

		// community : i escaped from community
		LUT_commu_esc_proba[H][c] = exp(-susc_commu_c * commu_transm_rate);

		// person-to-person : [1st] escaped from [2nd]
		for (int r = 0; r < c; ++r)
		{
			LUT_p2p_esc_proba[H][r][c] = exp(-house.group_rate(order[r], "Infectivity") * susc_house_c * p2p_transm_rate);
			LUT_p2p_esc_proba[H][c][r] = exp(-infect_c * house.group_rate(order[r], "Susc-House") * p2p_transm_rate);
		}

		// note : infector always escape himself (not relevant)
		LUT_p2p_esc_proba[H][c][c] = 1;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : initialize model  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the Lookup Table. */
void Model::init_model_for_inference()
{
	//==== Unknown status are not modeled =====//
	if (report_data("n.miss.stat") > 0)
	{
		tools::error("A missing status was found : incompatible with the model.");
	}

	//===== Allocate memory for the escape matrix =====//
	for (size_t H = 0; H < n_houses(); H++)
	{
		int h_size = house_size(H);
		LUT_commu_esc_proba.emplace_back(h_size, 1.0);
		LUT_p2p_esc_proba.emplace_back(h_size, std::vector<double>(h_size, 1.0));
	}

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_escape_patterns_lookuptable(report_data("max.n.HH.cases"));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////  Model : compute single house loglike  /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood for a single household. */
double Model::compute_house_loglike(const Household &house) const
{
	size_t H = house.index();
	int n_infec = house.n_infectious();

	//===== Fill Escape Probabilities =====//
	fill_house_escape_matrix(house);

	//=====compute the loglikelihood=====//
	return bax::chain_binomial_heter_loglikelihood(n_infec,
												   LUT_commu_esc_proba[H], 
												   LUT_p2p_esc_proba[H]);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// Model : Methods (simulation)  ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Simulate the infection within a household. */
void Model::simulate_single_house(Household &house)
{
	//=========================//
	//===== ESCAPE MATRIX =====//
	//=========================//

	size_t H = house.index();

	fill_house_escape_matrix(house);

	//==================================================//
	//===== GENERATION 0 : infection by community ======//
	//==================================================//

	for (int i = 0; i < house.size(); i++)
	{
		if (uniform() < LUT_commu_esc_proba[H][i])
		{
			house.change_status(i, SUSCEPTIBLE);
		}
		else
		{
			house.change_status(i, INFECTIOUS);
		}
	}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	std::set<int> old_gen_set = {};
	std::set<int> new_gen_set = house.infectious_set();

	while ((house.n_susceptible() > 0) && (!new_gen_set.empty()))
	{
		old_gen_set = new_gen_set;
		new_gen_set.clear();

		for (int i = 0; i < house.size(); i++)
		{
			if (house.status(i) == SUSCEPTIBLE)
			{
				for (int k : old_gen_set)
				{
					if (uniform() >= LUT_p2p_esc_proba[H][k][i])
					{
						house.change_status(i, INFECTIOUS);
						new_gen_set.insert(i);
					}
				}
			}
		}
	}
}
