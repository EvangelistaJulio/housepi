// model.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"
#include "binomial.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  GLOBAL VARIABLES  ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

//----- Attributes for Look-Up-Tables -----//
std::vector<std::vector<double>> LUT_susc_rates;
std::vector<std::vector<double>> LUT_infec_rates;
std::vector<std::vector<double>> LUT_commu_esc_proba;
std::vector<squarematrix> LUT_p2p_esc_proba;

//-----------------------------------------------------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

/* Assign the parameters controlling infectivity for each individual */
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

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* The person-to-person escape probabilities stored inside a matrix.
	Important: it must be re-ordered to have seropositives first. */
void Model::fill_house_escape_matrix(const Household &house) const
{
	//----- ordered members group rates : 1) infectious, 2) susceptibles -----//
	size_t H = house.index();
	int hh_size = house.size();
	int count = 0;

	double p2p_transm_rate = param("beta");
	double commu_transm_rate = param("alpha");

	// 1) ordered members group rates : infectious
	for (int i : house.infectious_set())
	{
		LUT_susc_rates[H][count] = house.group_rate(i, "Susceptibility");
		LUT_infec_rates[H][count] = house.group_rate(i, "Infectivity");
		++count;
	}

	// 2) ordered members group rates : susceptibles
	for (int i : house.susceptible_set())
	{
		LUT_susc_rates[H][count] = house.group_rate(i, "Susceptibility");
		LUT_infec_rates[H][count] = house.group_rate(i, "Infectivity");
		++count;
	}

	//----- assign values to matrix -----//
	for (int i = 0; i < hh_size; ++i)
	{
		// community : i escaped from community
		LUT_commu_esc_proba[H][i] = exp(-LUT_susc_rates[H][i] * commu_transm_rate);

		// person-to-person : i escaped from k
		for (int k = 0; k < i; ++k)
		{
			LUT_p2p_esc_proba[H][k][i] = exp(-LUT_susc_rates[H][i] * LUT_infec_rates[H][k] * p2p_transm_rate);
			LUT_p2p_esc_proba[H][i][k] = exp(-LUT_susc_rates[H][k] * LUT_infec_rates[H][i] * p2p_transm_rate);
		}

		// note : infector always escape himself (not relevant)
		LUT_p2p_esc_proba[H][i][i] = 1;
	}
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
	for (size_t H = 0; H < num_houses(); H++)
	{
		int h_size = house_size(H);
		LUT_susc_rates.emplace_back(h_size, 1.0);
		LUT_infec_rates.emplace_back(h_size, 1.0);
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
	int num_infec = house.num_infectious();

	//===== Fill Escape Probabilities =====//
	fill_house_escape_matrix(house);

	//=====compute the loglikelihood=====//
	return bax::chain_binomial_heter_loglikelihood(num_infec,
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
	int hh_size = house.size();

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

	while ((house.num_susceptible() > 0) && (!new_gen_set.empty()))
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
