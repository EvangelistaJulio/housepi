// model.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//-----------------------------------------------------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

/* Assign the parameters controlling infectivity for each individual */
void Model::assign_rate_per_group_to_individual(Individual &this_indiv)
{
	//=====Susceptibility=====//
	if (this_indiv.get_data("age") <= get_value("child_age"))
	{ // child
		this_indiv.assign_rate_per_group("Susceptibility", "sus_child");
	}
	else
	{ // adult
		if (this_indiv.get_data("sex") == 0)
		{ // adult women
			this_indiv.assign_rate_per_group("Susceptibility", "sus_women");
		}
		else
		{ // adult men : Susceptibility [REFERENCE]
			this_indiv.assign_rate_per_group("Susceptibility", REFERENCE);
		}
	}

	//=====Infectivity=====//
	if (this_indiv.get_data("age") > get_value("child_age"))
	{ // adult : Susceptiblity [REFERENCE]
		this_indiv.assign_rate_per_group("Infectivity", REFERENCE);
	}
	else
	{ // child
		this_indiv.assign_rate_per_group("Infectivity", "inf_child");
	}
}

/* The person-to-person escape probabilities stored inside a matrix.
	Important: it must be re-ordered to have seropositives first. */
void Model::fill_house_escape_matrix(int h,
									 std::vector<double> &commu_esc_proba,
									 std::vector<std::vector<double>> &p2p_esc_proba)
{
	//----- ordered members group rates : 1) infectious, 2) susceptibles -----//
	int hh_size = househ[h].house_size;
	std::vector<double> susc_rate(hh_size, 1.);
	std::vector<double> infec_rate(hh_size, 1.);
	int count = 0;

	// 1) ordered members group rates : infectious
	for (int i : househ[h].infectious_set)
	{
		susc_rate[count] = househ[h].get_group_rate(i, "Susceptibility");
		infec_rate[count] = househ[h].get_group_rate(i, "Infectivity");
		count++;
	}

	// 2) ordered members group rates : susceptibles
	for (int i : househ[h].susceptible_set)
	{
		susc_rate[count] = househ[h].get_group_rate(i, "Susceptibility");
		infec_rate[count] = househ[h].get_group_rate(i, "Infectivity");
		count++;
	}

	//----- assign values to matrix -----//
	for (int i = 0; i < hh_size; i++)
	{
		// community : i escaped from community
		commu_esc_proba[i] = exp(-susc_rate[i] * get_value("alpha"));

		// person-to-person : i escaped from k
		for (int k = 0; k < hh_size; k++)
		{
			if (k != i)
			{
				p2p_esc_proba[k][i] = exp(-susc_rate[i] * infec_rate[k] * get_value("beta"));
			}
			// else is equal to 1 (declared outside, not relevant)
		}
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
	// Unknown status are not modeled
	if (report_data["n.unk.stat"] > 0)
	{
		screen::raise_error("A missing status was found : incompatible with the model.");
	}

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_escape_patterns_lookuptable(report_data["max.n.HH.cases"]);
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////  Model : compute single house loglike  /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood for a single household. */
double Model::compute_house_loglike(int h)
{
	int hh_size = househ[h].house_size;
	int num_infec = househ[h].num_infectious;

	//=====Escape Matrix=====//
	std::vector<double> commu_esc_proba(hh_size, 1.);
	std::vector<std::vector<double>> p2p_esc_proba(hh_size, std::vector<double>(hh_size, 1.));
	fill_house_escape_matrix(h, commu_esc_proba, p2p_esc_proba);

	//=====compute the loglikelihood=====//
	return bax::chain_binomial_heter_loglikelihood(num_infec, commu_esc_proba, p2p_esc_proba);
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
void Model::simulate_single_house(int h)
{
	//=========================//
	//===== ESCAPE MATRIX =====//
	//=========================//

	int hh_size = househ[h].house_size;
	std::vector<double> commu_esc_proba(hh_size, 1.);
	std::vector<std::vector<double>> p2p_esc_proba(hh_size, std::vector<double>(hh_size, 1.));

	fill_house_escape_matrix(h, commu_esc_proba, p2p_esc_proba);

	//==================================================//
	//===== GENERATION 0 : infection by community ======//
	//==================================================//

	for (int i = 0; i < househ[h].house_size; i++)
	{
		if (uniform() < commu_esc_proba[i])
		{
			househ[h].change_status(i, SUSCEPTIBLE);
		}
		else
		{
			househ[h].change_status(i, INFECTIOUS);
		}
	}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	std::set<int> old_gen_set = {};
	std::set<int> new_gen_set = househ[h].infectious_set;

	while ((househ[h].num_susceptible > 0) && (!new_gen_set.empty()))
	{
		old_gen_set = new_gen_set;
		new_gen_set.clear();

		for (int i = 0; i < househ[h].house_size; i++)
		{
			if (househ[h].is_susceptible[i])
			{
				for (int k : old_gen_set)
				{
					if (uniform() >= p2p_esc_proba[k][i])
					{
						househ[h].change_status(i, INFECTIOUS);
						new_gen_set.insert(i);
					}
				}
			}
		}
	}
}
