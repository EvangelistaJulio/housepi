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

/* The house escape matrix is a flattern vector containing escape probability. */
void Model::fill_house_escape_matrix(int h, std::vector<double> &escape_matrix)
{
	int hh_size = househ[h].house_size;
	double roi;

	// ordered members group rates : 1) infectious, 2) susceptibles
	int count = 0;
	std::vector<double> susc_rates(1., hh_size), infec_rates(1., hh_size);
	for (int i : househ[h].infectious_list)
	{
		susc_rates[count] = househ[h].indiv[i].get_group_rate("Susceptibility");
		infec_rates[count] = househ[h].indiv[i].get_group_rate("Infectivity");
		count++;
	}
	for (int i : househ[h].susceptibles_list)
	{
		susc_rates[count] = househ[h].indiv[i].get_group_rate("Susceptibility");
		infec_rates[count] = househ[h].indiv[i].get_group_rate("Infectivity");
		count++;
	}

	// assign values to matrix
	for (int i = 0; i < hh_size; i++)
	{
		// community terms
		escape_matrix[i] = exp(-susc_rates[i] * get_value("alpha"));

		// person-to-person terms
		for (int k = 0; k < hh_size; k++)
		{
			if (k != i)
			{
				roi = get_value("beta") * susc_rates[i] * infec_rates[k];
				escape_matrix[i + (k + 1) * hh_size] = exp(-roi);
			}
			// else is equal to 1 (by declaration outside)
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
		raise_error_message("A missing status was found : incompatible with the model.");
	}

	//=====Initialize look-up-table to largest possible value=====//
	int max_num = report_data["max.n.HH.cases"];
	bax::init_combinations_lookuptable(max_num);
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
	std::vector<double> escape_matrix((hh_size + 1) * hh_size, 1.);
	fill_house_escape_matrix(h, escape_matrix);

	//=====compute the loglikelihood=====//
	// initialize array of indexes to {1, ..., hh_size}
	bax::vec_ushort indidx(hh_size);
	for (ushort i = 1; i <= hh_size; i++)
	{
		indidx[i - 1] = i;
	}

	double likelihood = bax::ll_chain_binomial_heterogeneous(
		hh_size, num_infec, indidx, escape_matrix);

	return log(likelihood);
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
	std::vector<double> escape_matrix((hh_size + 1) * hh_size, 1.);
	fill_house_escape_matrix(h, escape_matrix);

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
		if (uniform() > escape_matrix[i])
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
					if (uniform() > escape_matrix[i + (k + 1) * hh_size])
					{
						househ[h].change_status(i, INFECTIOUS);
						new_infectors_list.push_back(i);
					}
				}
			}
		}
	}
}
