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
	if (this_indiv.get_data("age") <= get_param("child_age"))
	{ // child
		this_indiv.assign_rate_per_group("Susceptibility", "s_child");
	}
	else
	{ // adult
		if (this_indiv.get_data("sex") == 0)
		{ // adult women
			this_indiv.assign_rate_per_group("Susceptibility", "s_women");
		}
		else
		{ // adult men : Susceptibility [REFERENCE]
			this_indiv.assign_rate_per_group("Susceptibility", REFERENCE);
		}
	}

	//=====Infectivity=====//
	if (this_indiv.get_data("age") > get_param("child_age"))
	{ // adult : Susceptiblity [REFERENCE]
		this_indiv.assign_rate_per_group("Infectivity", REFERENCE);
	}
	else
	{ // child
		this_indiv.assign_rate_per_group("Infectivity", "i_child");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_rate_per_contact_to_household(int h, int i, int k)
{
	// from `i` -> to `k`
	// children if younger than 12 y.o., adults otherwise
	bool i_is_child = (get_data(h, i, "age") < get_param("child_age"));
	bool i_is_adult = (i_is_child == false);

	bool k_is_child = (get_data(h, k, "age") < get_param("child_age"));
	bool k_is_adult = (k_is_child == false);

	if (i_is_child && k_is_child)
	{ // Child to Child
		house[h].assign_rate_per_contact(i, k, "c_kid_to_kid");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		house[h].assign_rate_per_contact(i, k, "c_kid_to_adu");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		house[h].assign_rate_per_contact(i, k, "c_adu_to_kid");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		house[h].assign_rate_per_contact(i, k, REFERENCE);
	}
}

/* The person-to-person escape probabilities stored inside a matrix.
	Important: it must be re-ordered to have seropositives first. */
void Model::fill_house_escape_matrix(int h, vector<double> &commu_esc_proba, squarematrix &p2p_esc_proba)
{
	//----- ordered members group rates : 1) infectious, 2) susceptibles -----//
	int hh_size = house[h].house_size;
	vector<double> susc_rate(hh_size, 1.);
	vector<double> infec_rate(hh_size, 1.);
	int count = 0;

	// 1) ordered members group rates : infectious
	for (int i : house[h].infectious_set)
	{
		susc_rate[count] = house[h].get_group_rate(i, "Susceptibility");
		infec_rate[count] = house[h].get_group_rate(i, "Infectivity");
		count++;
	}

	// 2) ordered members group rates : susceptibles
	for (int i : house[h].susceptible_set)
	{
		susc_rate[count] = house[h].get_group_rate(i, "Susceptibility");
		infec_rate[count] = house[h].get_group_rate(i, "Infectivity");
		count++;
	}

	//----- assign values to matrix -----//
	for (int i = 0; i < hh_size; i++)
	{
		// community : i escaped from community
		commu_esc_proba[i] = exp(-susc_rate[i] * get_param("alpha"));

		// person-to-person : i escaped from k
		for (int k = 0; k < hh_size; k++)
		{
			if (k != i)
			{
				p2p_esc_proba[k][i] = exp(-susc_rate[i] * get_contact_rate(h, k, i) * infec_rate[k] * get_param("beta"));
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
		basics::raise_error("A missing status was found : incompatible with the model.");
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
	int hh_size = house[h].house_size;
	int num_infec = house[h].num_infectious;

	//=====Escape Matrix=====//
	vector<double> commu_esc_proba(hh_size, 1.);
	vector<vector<double>> p2p_esc_proba(hh_size, vector<double>(hh_size, 1.));
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

	int hh_size = house[h].house_size;
	vector<double> commu_esc_proba(hh_size, 1.);
	vector<vector<double>> p2p_esc_proba(hh_size, vector<double>(hh_size, 1.));

	fill_house_escape_matrix(h, commu_esc_proba, p2p_esc_proba);

	//==================================================//
	//===== GENERATION 0 : infection by community ======//
	//==================================================//

	for (int i = 0; i < house[h].house_size; i++)
	{
		if (uniform() < commu_esc_proba[i])
		{
			house[h].change_status(i, SUSCEPTIBLE);
		}
		else
		{
			house[h].change_status(i, INFECTIOUS);
		}
	}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	std::set<int> old_gen_set = {};
	std::set<int> new_gen_set = house[h].infectious_set;

	while ((house[h].num_susceptible > 0) && (!new_gen_set.empty()))
	{
		old_gen_set = new_gen_set;
		new_gen_set.clear();

		for (int i = 0; i < house[h].house_size; i++)
		{
			if (house[h].is_susceptible[i])
			{
				for (int k : old_gen_set)
				{
					if (uniform() >= p2p_esc_proba[k][i])
					{
						house[h].change_status(i, INFECTIOUS);
						new_gen_set.insert(i);
					}
				}
			}
		}
	}
}
