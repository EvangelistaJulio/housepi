// model.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  * *** * CHAIN BINOMIAL HETEROGENEOUS * *** *  /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//------------------------------------------------------------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- \/ -------------------------------------------------------//

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
	Important: it must be re-ordered to have seroposatives first. */
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
	for (int i : househ[h].infectious_list)
	{
		susc_rate[count] = househ[h].indiv[i].get_group_rate("Susceptibility");
		infec_rate[count] = househ[h].indiv[i].get_group_rate("Infectivity");
		count++;
	}

	// 2) ordered members group rates : susceptibles
	for (int i : househ[h].susceptibles_list)
	{
		susc_rate[count] = househ[h].indiv[i].get_group_rate("Susceptibility");
		infec_rate[count] = househ[h].indiv[i].get_group_rate("Infectivity");
		count++;
	}

	//----- assign values to matrix -----//
	double aux;
	for (int i = 0; i < hh_size; i++)
	{
		// community : i escaped from community
		commu_esc_proba[i] = exp(-susc_rate[i] * get_value("alpha"));

		// person-to-person : i escaped from k
		for (int k = 0; k < hh_size; k++)
		{
			if (k != i)
			{
				aux = get_value("beta") * std::pow(get_value("typic_size") / hh_size, get_value("gamma"));
				p2p_esc_proba[k][i] = exp(-susc_rate[i] * infec_rate[k] * aux);
			}
			// else is equal to 1 (declared outside, not relevant)
		}
	}
}

//------------------------------------------------- /\ -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//------------------------------------------------- || -------------------------------------------------------//
//---------------------------------------- - - - CHANGE-ME - - - ---------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////  Model : initialize model  //////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the Lookup Table. */
void Model::initialize_model()
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
double Model::compute_single_house_loglike(int h)
{
	int hh_size = househ[h].house_size;
	int num_infec = househ[h].infectious_list.size();

	//=====Escape Matrix=====//
	std::vector<double> commu_esc_proba(hh_size, 1.);
	std::vector<std::vector<double>> p2p_esc_proba(hh_size, std::vector<double>(hh_size, 1.));
	fill_house_escape_matrix(h, commu_esc_proba, p2p_esc_proba);

	//=====compute the loglikelihood=====//
	return bax::chain_binomial_heter_loglikelihood(num_infec, commu_esc_proba, p2p_esc_proba);
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (simulation)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Simulate the infection within a household. */
void Model::simulate_single_house(int h)
{
	int hh_size = househ[h].house_size;
	std::vector<int> store_old_list = {};
	std::vector<int> old_infectors_list = {}, new_infectors_list = {};

	//=====Escape Matrix=====//
	std::vector<double> commu_esc_proba(hh_size, 1.);
	std::vector<std::vector<double>> p2p_esc_proba(hh_size, std::vector<double>(hh_size, 1.));
	fill_house_escape_matrix(h, commu_esc_proba, p2p_esc_proba);

	//=====INITIALIZE=====//
	// set all cases to healthy (change infectious to susceptibles)
	store_old_list = househ[h].infectious_list;
	for (int i : store_old_list)
	{
		househ[h].change_status(i, SUSCEPTIBLE);
	}

	//====GENERATION 0 : infection by community=====//
	store_old_list = househ[h].susceptibles_list;
	for (int i : store_old_list)
	{
		if (uniform() > commu_esc_proba[i])
		{
			househ[h].change_status(i, INFECTIOUS);
			new_infectors_list.push_back(i);
		}
	}

	//====GENERATION >=1 : infection by previous generatioon=====//
	while ((househ[h].susceptibles_list.size() > 0) && (new_infectors_list.size() > 0))
	{
		old_infectors_list = new_infectors_list;
		new_infectors_list.clear();

		store_old_list = househ[h].susceptibles_list;
		for (int i : store_old_list)
		{
			for (int k : old_infectors_list)
			{
				if (househ[h].indiv[i].get_status() == SUSCEPTIBLE)
				{
					if (uniform() > p2p_esc_proba[k][i])
					{
						househ[h].change_status(i, INFECTIOUS);
						new_infectors_list.push_back(i);
					}
				}
			}
		}
	}
}
