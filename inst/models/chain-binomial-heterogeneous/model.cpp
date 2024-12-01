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

//-----------------------------------------------------------------------------------------------//
//---------------------------------------- CHANGE-ME --------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

/* [CHANGE-ME] Assign the parameters controlling infectivity for each individual. */
void Model::assign_group_rate_to_individual(Household &house, int i)
{
	std::string rateID;
	bool is_adult = house.data(i, "age") > param("child_age");
	bool is_male = (house.data(i, "sex") == 1);

	//===== Infectivity in Community =====//
	rateID = "Infectivity";
	if (is_adult)
		house.assign_group_rate(i, rateID, REFERENCE);
	else
		house.assign_group_rate(i, rateID, "inf_child");

	//===== Susceptibility in Community =====//
	rateID = "Susc-Community";
	if (is_adult)
		if (is_male)
			house.assign_group_rate(i, rateID, REFERENCE);
		else
			house.assign_group_rate(i, rateID, "comsus_woman");
	else
		house.assign_group_rate(i, rateID, "comsus_child");

	//===== Susceptibility in House =====//
	rateID = "Susc-House";
	if (is_adult)
		if (is_male)
			house.assign_group_rate(i, rateID, REFERENCE);
		else
			house.assign_group_rate(i, rateID, "sus_woman");
	else
		house.assign_group_rate(i, rateID, "sus_child");
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
////////////////////////////////////  Model : initialize model  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the Lookup Table. */
void Model::init_model_for_simulation()
{
	//===== Allocate memory for the escape matrix =====//
	assign_extra_vector_to_all_houses({"commu-escprob"});
	assign_extra_matrix_to_all_houses({"p2p-escprob"});
}

/* Initialize the Lookup Table. */
void Model::init_model_for_inference()
{
	//==== Unknown status are not modeled =====//
	if (report_data("n.miss.stat") > 0)
	{
		tools::error("A missing status was found : incompatible with the model.");
	}

	//===== Allocate memory for the escape matrix =====//
	init_model_for_simulation();

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_escape_patterns_lookuptable(report_data("max.n.HH.cases"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Model : fill house escape matrix  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The person-to-person escape probabilities stored inside a matrix.
	Important: it must be re-ordered to have seropositives first. */
void Model::fill_house_escape_matrix(Household &house) const
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
	for (int c = 0; c < house.size(); ++c)
	{
		double susc_commu_c = house.group_rate(order[c], "Susc-Community");
		double infect_c = house.group_rate(order[c], "Infectivity");
		double susc_house_c = house.group_rate(order[c], "Susc-House");

		// community : i escaped from community
		house.set_extra_vector("commu-escprob", c, exp(-susc_commu_c * commu_transm_rate));

		// person-to-person : [1st] escaped from [2nd]
		for (int r = 0; r < c; ++r)
		{
			house.set_extra_matrix("p2p-escprob", r, c,
								   exp(-house.group_rate(order[r], "Infectivity") * susc_house_c * p2p_transm_rate));
			house.set_extra_matrix("p2p-escprob", c, r,
								   exp(-infect_c * house.group_rate(order[r], "Susc-House") * p2p_transm_rate));
		}

		// note : infector always escape himself (not relevant)
		house.set_extra_matrix("p2p-escprob", c, c, 1.0);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////  Model : compute single house loglike  /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood for a single household. */
double Model::compute_house_loglike(Household &house) const
{
	//===== Fill Escape Probabilities =====//
	fill_house_escape_matrix(house);

	//=====compute the loglikelihood=====//
	return bax::chain_binomial_heter_loglikelihood(house.n_infectious(),
												   house.extra_vector("commu-escprob"),
												   house.extra_matrix("p2p-escprob"));
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

	fill_house_escape_matrix(house);

	//==================================================//
	//===== GENERATION 0 : infection by community ======//
	//==================================================//

	for (int i = 0; i < house.size(); i++)
	{
		if (uniform() < house.extra_vector("commu-escprob", i))
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
					if (uniform() >= house.extra_matrix("commu-escprob", k, i))
					{
						house.change_status(i, INFECTIOUS);
						new_gen_set.insert(i);
					}
				}
			}
		}
	}
}
