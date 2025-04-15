/*
 * model.cpp
 *
 * Model template of the 'housepi' package.
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
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

/* Assign all group rates to each individual. */
void Model::assign_group_rate_to_individual(Household &house, indiv i)
{
	rateID r;
	UNUSED(house);
	UNUSED(i);

	//----- community -----//
	r = "community";
	emplace_group_rate(r, REFERENCE);

	//----- infectivity -----//
	r = "infectivity";
	emplace_group_rate(r, REFERENCE);

	//----- susceptibility -----//
	r = "susceptibility";
	emplace_group_rate(r, REFERENCE);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : initialize model  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the Lookup Table. */
void Model::init_model_for_inference()
{
	//==== Unknown status are not modeled =====//
	if (get_report_indiv("missing") > 0)
		tools::error(
			"A missing status was found : incompatible with the model.");

	//===== Allocate memory for the escape matrix =====//
	allocate_extra_house_integer_vector("order");
	allocate_extra_house_double_vector("commu-escprob");
	allocate_extra_house_double_matrix("p2p-escprob");

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_escape_patterns_lookuptable(get_report_hh_size("max cases").first);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  Model : fill house escape matrix  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The person-to-person escape probabilities stored inside a matrix.
	Important: it must be re-ordered to have seropositives first. */
void Model::fill_house_escape_matrix(Household &house) const
{
	const double commu_transm_rate = community_risk();
	const double p2p_transm_rate = person_to_person_risk();

	//----- ordered members group rates : 1) infectious, 2) susceptibles -----//
	auto &order = house.extra_integer_vector("order");

	indiv m = 0;
	for (indiv i : house.infected_status_list())
		order[m++] = i;
	for (indiv s : house.healthy_status_list())
		order[m++] = s;

	//----- Assign values to LUT matrix -----//
	auto &commu_escprob = house.extra_double_vector("commu-escprob");
	auto &p2p_escprob = house.extra_double_matrix("p2p-escprob");

	for (indiv c = 0; c < house.size(); ++c)
	{
		indiv idx_c = order[c];
		double infect_c = house.group_rate(idx_c, "infectivity");
		double susc_house_c = house.group_rate(idx_c, "susceptibility");

		//----- community -----//
		commu_escprob[c] = exp(-commu_transm_rate * house.group_rate(idx_c, "community"));

		//----- person-to-person -----//
		for (indiv r = 0; r < c; ++r)
		{
			indiv idx_r = order[r];

			// [c] : infector -> [r] : infectee
			p2p_escprob[c][r] = exp(-infect_c * p2p_transm_rate * house.group_rate(idx_r, "susceptibility"));

			// [r] : infector -> [c] : infectee
			p2p_escprob[r][c] = exp(-house.group_rate(idx_r, "infectivity") * p2p_transm_rate * susc_house_c);
		}

		// note : infector always "escapes" himself (not relevant)
		p2p_escprob[c][c] = 1.0;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////  Model : compute single house loglike  /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood for a single household. */
double Model::compute_house_obs_loglike(Household &house) const
{
	//===== Fill Escape Probabilities =====//
	fill_house_escape_matrix(house);

	//=====compute the loglikelihood=====//
	const double loglikelihood = bax::chain_binomial_heter_loglikelihood(
		house.n_infec_status(),
		house.extra_double_vector("commu-escprob"),
		house.extra_double_matrix("p2p-escprob")
	);
	return loglikelihood * house.multiplicity();
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
// UPDATE : collapse to fill_house_escape_matrix(house)
void Model::simulate_single_house(Household &house)
{
	const double commu_transm_rate = community_risk();
	const double p2p_transm_rate = person_to_person_risk();

	//=================================================//
	//===== GENERATION 0 : infection by community =====//
	//=================================================//

	std::vector<indiv> new_gen_set = {};

	for (indiv i = 0; i < house.size(); ++i)
		if (uniform() < exp(-commu_transm_rate * house.group_rate(i, "community")))
			house.change_status(i, HEALTHY_STATUS);
		else
		{
			house.change_status(i, INFECTED_STATUS);
			new_gen_set.push_back(i);
		}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	std::vector<indiv> old_gen_set = {};

	while ((house.n_healthy_status() > 0) && (!new_gen_set.empty()))
	{
		std::swap(old_gen_set, new_gen_set);
		new_gen_set.clear();

		// note : loop over 'house.healthy_status_list()' is not allowed because of 'change_status'
		for (indiv i = 0; i < house.size(); ++i)
			if (house.status(i) == HEALTHY_STATUS)
			{
				double sus_i = house.group_rate(i, "susceptibility");

				double r = 0.;
				for (indiv k : old_gen_set)
					r += house.group_rate(k, "infectivity");

				if (uniform() >= exp(-r * p2p_transm_rate * sus_i))
				{
					house.change_status(i, INFECTED_STATUS);
					new_gen_set.push_back(i);
				}
			}
	}
}
