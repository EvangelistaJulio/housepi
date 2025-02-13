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

/* Assign all group rates for each individual. */
void Model::assign_group_rate_to_individual(Household &house, indiv i)
{
	rateID r;
	bool is_adult = house.data(i, "age") > param("child_age");
	bool is_male = (house.data(i, "sex") == 1);

	//----- infectivity -----//
	r = "infectivity";
	if (is_adult)
		emplace_group_rate(r, REFERENCE);
	else
		emplace_group_rate(r, "inf_child");

	//----- community_susceptibility -----//
	r = "community_susceptibility";
	if (is_adult && is_male)
		emplace_group_rate(r, REFERENCE);
	else if (is_adult && (!is_male))
		emplace_group_rate(r, "com_sus_woman");
	else
		emplace_group_rate(r, "com_sus_child");

	//----- susceptibility -----//
	r = "susceptibility";
	if (is_adult && is_male)
		emplace_group_rate(r, REFERENCE);
	else if (is_adult && (!is_male))
		emplace_group_rate(r, "sus_woman");
	else
		emplace_group_rate(r, "sus_child");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Model : initialize model  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the Lookup Table. */
void Model::init_model_for_inference()
{
	//==== Unknown status are not modeled =====//
	if (report_data("n.miss.stat") > 0)
		tools::error(
			"A missing status was found : incompatible with the model.");

	//===== Allocate memory for the escape matrix =====//
	assign_extra_house_integer_vector("order");
	assign_extra_house_double_vector("commu-escprob");
	assign_extra_house_double_matrix("p2p-escprob");

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
	const double commu_transm_rate = community_risk();
	const double p2p_transm_rate = person_to_person_risk();

	//----- ordered members group rates : 1) infectious, 2) susceptibles -----//
	auto &order = house.extra_integer_vector("order");

	indiv m = 0;
	for (indiv i : house.infectious_set())
		order[m++] = i;
	for (indiv s : house.susceptible_set())
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
		commu_escprob[c] = exp(-commu_transm_rate * house.group_rate(idx_c, "community_susceptibility"));

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
double Model::compute_house_loglike(Household &house) const
{
	//===== Fill Escape Probabilities =====//
	fill_house_escape_matrix(house);

	//=====compute the loglikelihood=====//
	return bax::chain_binomial_heter_loglikelihood(house.n_infectious(),
												   house.extra_double_vector("commu-escprob"),
												   house.extra_double_matrix("p2p-escprob"));
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

	//=========================//
	//===== ESCAPE MATRIX =====//
	//=========================//

	for (indiv i = 0; i < house.size(); ++i)
		if (uniform() < exp(-commu_transm_rate * house.group_rate(i, "community_susceptibility")))
			house.change_status(i, HEALTHY_STATUS);
		else
			house.change_status(i, INFECTED_STATUS);

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	std::set<indiv> old_gen_set = {};
	std::set<indiv> new_gen_set = house.infectious_set();

	while ((house.n_susceptible() > 0) && (!new_gen_set.empty()))
	{
		std::swap(old_gen_set, new_gen_set);
		new_gen_set.clear();

		// note : loop over 'house.susceptible_set()' is not allowed because of 'change_status'
		for (indiv i = 0; i < house.size(); ++i)
		{
			if (house.status(i) == HEALTHY_STATUS)
			{
				double sus_i = house.group_rate(i, "susceptibility");

				double r = 0.;
				for (indiv k : old_gen_set)
					r += house.group_rate(k, "infectivity");

				if (uniform() >= exp(-r * p2p_transm_rate * sus_i))
				{
					house.change_status(i, INFECTED_STATUS);
					new_gen_set.insert(i);
				}
			}
		}
	}
}
