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
void Model::assign_group_rate_to_individual(Individual &member)
{
	rateID r;
	UNUSED(member);

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
	allocate_cache_house_integer_vector("order");
	allocate_cache_house_double_vector("commu-escprob");
	allocate_cache_house_double_matrix("p2p-escprob");

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
	auto &order = house.cache_integer_vector("order");

	indiv m = 0;
	for (indiv i : house.get_infected_status_hhindexes())
		order[m++] = i;
	for (indiv s : house.get_healthy_status_hhindexes())
		order[m++] = s;

	//----- Assign values to LUT matrix -----//
	auto &commu_escprob = house.cache_double_vector("commu-escprob");
	auto &p2p_escprob = house.cache_double_matrix("p2p-escprob");

	for (indiv c = 0; c < house.size(); ++c)
	{
		Individual &member_c = house.get_member(order[c]);
		double infect_c = member_c.group_rate("infectivity");
		double susc_house_c = member_c.group_rate("susceptibility");

		//----- community -----//
		commu_escprob[c] = exp(-commu_transm_rate * member_c.group_rate("community"));

		//----- person-to-person -----//
		for (indiv r = 0; r < c; ++r)
		{
			Individual &member_r = house.get_member(order[r]);

			// [c] : infector -> [r] : infectee
			p2p_escprob[c][r] = exp(-infect_c * p2p_transm_rate * member_r.group_rate("susceptibility"));

			// [r] : infector -> [c] : infectee
			p2p_escprob[r][c] = exp(-member_r.group_rate("infectivity") * p2p_transm_rate * susc_house_c);
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
		house.cache_double_vector("commu-escprob"),
		house.cache_double_matrix("p2p-escprob"));
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

	std::vector<Individual*> new_gen_set = {};

	for (Individual &member : house)
		if (uniform() < exp(-commu_transm_rate * member.group_rate("community")))
			house.change_status(member, HEALTHY_STATUS);
		else
		{
			house.change_status(member, INFECTED_STATUS);
			new_gen_set.push_back(&member);
		}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	std::vector<Individual*> old_gen_set = {};

	while ((house.n_healthy_status() > 0) && (!new_gen_set.empty()))
	{
		std::swap(old_gen_set, new_gen_set);
		new_gen_set.clear();

		for (Individual &member_i : house.healthy())
		{
			double sus_i = member_i.group_rate("susceptibility");

			double r = 0.;
			for (Individual *member_k : old_gen_set)
				r += member_k->group_rate("infectivity");

			if (uniform() >= exp(-r * p2p_transm_rate * sus_i))
			{
				house.change_status(member_i, INFECTED_STATUS);
				new_gen_set.push_back(&member_i);
			}
		}
	}
}
