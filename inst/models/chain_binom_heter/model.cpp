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
void Model::assign_group_rate_to_individual(Household &this_house, int i)
{
	//=====Susceptibility=====//
	if (this_house.get_data(i, "age") <= get_param("child_age"))
	{ // child
		this_house.assign_group_rate(i, "Susceptibility", "sus_child");
	}
	else
	{ // adult
		if (this_house.get_data(i, "sex") == 0)
		{ // adult women
			this_house.assign_group_rate(i, "Susceptibility", "sus_women");
		}
		else
		{ // adult men : Susceptibility [REFERENCE]
			this_house.assign_group_rate(i, "Susceptibility", REFERENCE);
		}
	}

	//=====Infectivity=====//
	if (this_house.get_data(i, "age") > get_param("child_age"))
	{ // adult : Susceptiblity [REFERENCE]
		this_house.assign_group_rate(i, "Infectivity", REFERENCE);
	}
	else
	{ // child
		this_house.assign_group_rate(i, "Infectivity", "inf_child");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_contact_rate_to_household(Household &this_house, int from, int to)
{
	// from `i` -> to `k`
	// children if younger than 12 y.o., adults otherwise
	bool i_is_child = (this_house.get_data(from, "age") < get_param("child_age"));
	bool i_is_adult = (i_is_child == false);

	bool k_is_child = (this_house.get_data(to, "age") < get_param("child_age"));
	bool k_is_adult = (k_is_child == false);

	if (i_is_child && k_is_child)
	{ // Child to Child
		this_house.assign_contact_rate(from, to, "cont_child_child");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		this_house.assign_contact_rate(from, to, "cont_child_adult");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		this_house.assign_contact_rate(from, to, "cont_adult_child");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		this_house.assign_contact_rate(from, to, REFERENCE);
	}
}

/* The person-to-person escape probabilities stored inside a matrix.
	Important: it must be re-ordered to have seropositives first. */
void Model::fill_house_escape_matrix(const Household &this_house, std::vector<double> &commu_esc_proba, squarematrix &p2p_esc_proba) const
{
	//----- ordered members group rates : 1) infectious, 2) susceptibles -----//
	int hh_size = this_house.size();
	std::vector<double> susc_rate(hh_size, 1.);
	std::vector<double> infec_rate(hh_size, 1.);
	int count = 0;

	// 1) ordered members group rates : infectious
	for (int i : this_house.get_infectious_set())
	{
		susc_rate[count] = this_house.get_group_rate(i, "Susceptibility");
		infec_rate[count] = this_house.get_group_rate(i, "Infectivity");
		count++;
	}

	// 2) ordered members group rates : susceptibles
	for (int i : this_house.get_susceptible_set())
	{
		susc_rate[count] = this_house.get_group_rate(i, "Susceptibility");
		infec_rate[count] = this_house.get_group_rate(i, "Infectivity");
		count++;
	}

	//----- assign values to matrix -----//
	for (int i = 0; i < hh_size; ++i)
	{
		// community : i escaped from community
		commu_esc_proba[i] = exp(-susc_rate[i] * get_param("alpha"));

		// person-to-person : i escaped from k
		// note : else is equal to 1 (declared outside, not relevant)
		for (int k = 0; k < i; ++k)
		{
			p2p_esc_proba[k][i] = exp(-susc_rate[i] * this_house.get_contact_rate(k, i) * infec_rate[k] * get_param("beta"));
			p2p_esc_proba[i][k] = exp(-susc_rate[k] * this_house.get_contact_rate(i, k) * infec_rate[i] * get_param("beta"));
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
	if (report_data("n.miss.stat") > 0)
	{
		tools::error("A missing status was found : incompatible with the model.");
	}

	//=====Initialize look-up-table to largest possible value=====//
	bax::init_escape_patterns_lookuptable(report_data("max.n.HH.cases"));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////  Model : compute single house loglike  /////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood for a single household. */
double Model::compute_house_loglike(const Household &this_house) const
{
	int hh_size = this_house.size();
	int num_infec = this_house.get_num_infectious();

	//=====Escape Matrix=====//
	std::vector<double> commu_esc_proba(hh_size, 1.);
	squarematrix p2p_esc_proba(hh_size, std::vector<double>(hh_size, 1.));
	fill_house_escape_matrix(this_house, commu_esc_proba, p2p_esc_proba);

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
void Model::simulate_single_house(Household &this_house)
{
	//=========================//
	//===== ESCAPE MATRIX =====//
	//=========================//

	int hh_size = this_house.size();
	std::vector<double> commu_esc_proba(hh_size, 1.);
	squarematrix p2p_esc_proba(hh_size, std::vector<double>(hh_size, 1.));

	fill_house_escape_matrix(this_house, commu_esc_proba, p2p_esc_proba);

	//==================================================//
	//===== GENERATION 0 : infection by community ======//
	//==================================================//

	for (int i = 0; i < this_house.size(); i++)
	{
		if (uniform() < commu_esc_proba[i])
		{
			this_house.change_status(i, SUSCEPTIBLE);
		}
		else
		{
			this_house.change_status(i, INFECTIOUS);
		}
	}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	std::set<int> old_gen_set = {};
	std::set<int> new_gen_set = this_house.get_infectious_set();

	while ((this_house.get_num_susceptible() > 0) && (!new_gen_set.empty()))
	{
		old_gen_set = new_gen_set;
		new_gen_set.clear();

		for (int i = 0; i < this_house.size(); i++)
		{
			if (this_house.get_status(i) == SUSCEPTIBLE)
			{
				for (int k : old_gen_set)
				{
					if (uniform() >= p2p_esc_proba[k][i])
					{
						this_house.change_status(i, INFECTIOUS);
						new_gen_set.insert(i);
					}
				}
			}
		}
	}
}
