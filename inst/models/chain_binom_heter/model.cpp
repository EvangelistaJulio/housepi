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
void Model::fill_house_escape_matrix(const Household &this_house) const
{
	//----- ordered members group rates : 1) infectious, 2) susceptibles -----//
	size_t H = this_house.get_index();
	int hh_size = this_house.size();
	int count = 0;

	double p2p_transm_rate = get_param("beta");
	double commu_transm_rate = get_param("alpha");

	// 1) ordered members group rates : infectious
	for (int i : this_house.get_infectious_set())
	{
		LUT_susc_rates[H][count] = this_house.get_group_rate(i, "Susceptibility");
		LUT_infec_rates[H][count] = this_house.get_group_rate(i, "Infectivity");
		++count;
	}

	// 2) ordered members group rates : susceptibles
	for (int i : this_house.get_susceptible_set())
	{
		LUT_susc_rates[H][count] = this_house.get_group_rate(i, "Susceptibility");
		LUT_infec_rates[H][count] = this_house.get_group_rate(i, "Infectivity");
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
			LUT_p2p_esc_proba[H][k][i] = exp(-LUT_susc_rates[H][i] * this_house.get_contact_rate(k, i) * LUT_infec_rates[H][k] * p2p_transm_rate);
			LUT_p2p_esc_proba[H][i][k] = exp(-LUT_susc_rates[H][k] * this_house.get_contact_rate(i, k) * LUT_infec_rates[H][i] * p2p_transm_rate);
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
	for (size_t H = 0; H < get_num_house(); H++)
	{
		int h_size = house(H).size();
		LUT_susc_rates.emplace_back(h_size, 1);
		LUT_infec_rates.emplace_back(h_size, 1);
		LUT_commu_esc_proba.emplace_back(h_size, 1);
		LUT_p2p_esc_proba.emplace_back(h_size, std::vector<double>(h_size, 1.));
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
	size_t H = this_house.get_index();
	int num_infec = this_house.get_num_infectious();

	//===== Fill Escape Probabilities =====//
	fill_house_escape_matrix(this_house);

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
void Model::simulate_single_house(Household &this_house)
{
	//=========================//
	//===== ESCAPE MATRIX =====//
	//=========================//

	size_t H = this_house.get_index();
	int hh_size = this_house.size();

	fill_house_escape_matrix(this_house);

	//==================================================//
	//===== GENERATION 0 : infection by community ======//
	//==================================================//

	for (int i = 0; i < this_house.size(); i++)
	{
		if (uniform() < LUT_commu_esc_proba[H][i])
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
					if (uniform() >= LUT_p2p_esc_proba[H][k][i])
					{
						this_house.change_status(i, INFECTIOUS);
						new_gen_set.insert(i);
					}
				}
			}
		}
	}
}
