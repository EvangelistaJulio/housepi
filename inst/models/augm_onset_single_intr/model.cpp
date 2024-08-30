/*
 * model.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////  * *** * AUGMENTED ONSET - SINGLE INTRODUCTION * *** *  /////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "model.h"

//-----------------------------------------------------------------------------------------------//
//---------------------------------- - - - CHANGE-ME - - - --------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- || ------------------------------------------------//
//------------------------------------------- \/ ------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate per group to individual  ////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

void Model::assign_rate_per_group_to_individual(Individual &this_indiv)
{
	//===== rate n.0 : Infectivity =====//
	if (this_indiv.get_data("vacc_status") == 0)
	{ // naive : Infectivity [REFERENCE]
		this_indiv.assign_rate_per_group("Infectivity", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		this_indiv.assign_rate_per_group("Infectivity", "i_vacc");
	}

	//===== rate n.1 : Susceptibility =====//
	if (this_indiv.get_data("vacc_status") == 0)
	{ // naive : Susceptibility [REFERENCE]
		this_indiv.assign_rate_per_group("Susceptibility", REFERENCE);
	}
	else
	{ // fully vaccinated : Infectivity
		this_indiv.assign_rate_per_group("Susceptibility", "s_vacc");
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
	bool i_is_child = (househ[h].indiv[i].get_data("age") < get_value("child_age"));
	bool i_is_adult = (i_is_child == false);

	bool k_is_child = (househ[h].indiv[k].get_data("age") < get_value("child_age"));
	bool k_is_adult = (k_is_child == false);

	if (i_is_child && k_is_child)
	{ // Child to Child
		househ[h].assign_rate_per_contact(i, k, "c_kid_to_kid");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		househ[h].assign_rate_per_contact(i, k, "c_kid_to_adu");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		househ[h].assign_rate_per_contact(i, k, "c_adu_to_kid");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		househ[h].assign_rate_per_contact(i, k, REFERENCE);
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
///////////////////////////////////  Model : get infection time  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
double Model::get_infec_time(int h, int i) const
{
	return househ[h].indiv[i].get_augm("onset_time") - househ[h].indiv[i].get_augm("incub_period");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : instantaneous infectivity profile  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Instanteneous infectivity profile. */
double Model::compute_infec_prof_inst(int h, int i, double t_transm)
{
	// tost : time from onset of symptoms to transmission
	double tost = t_transm - househ[h].indiv[i].get_augm("onset_time");
	double incub_period = househ[h].indiv[i].get_augm("incub_period");
	return LUT::infec_prof_func.pdf(tost, incub_period);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// Model : cumulative infectivity profile  ///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Cumulative infectivity profile. */
double Model::compute_infec_prof_cumul(int h, int i, double t_transm)
{
	// tost : time from onset of symptoms to transmission
	double tost = t_transm - househ[h].indiv[i].get_augm("onset_time");
	double incub_period = househ[h].indiv[i].get_augm("incub_period");
	return LUT::infec_prof_func.cdf(tost, incub_period);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : store dependencies of augmented variables  ///////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Store cumulative and instantaneous infectivity profiles according to an infector. */
void Model::augm_var_store_dependecies(int h, int i)
{
	double t_infec_i = get_infec_time(h, i);
	double t_infec_k;

	//============================//
	//=====Infectious Members=====//
	//============================//

	for (int k : househ[h].infectious_set)
	{
		t_infec_k = get_infec_time(h, k);

		if (t_infec_i < t_infec_k)
		{ // `i` : infector, `k` : infectee
			LUT::inst_profile[h][i][k] = compute_infec_prof_inst(h, i, t_infec_k);
			LUT::cumul_profile[h][i][k] = compute_infec_prof_cumul(h, i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			LUT::inst_profile[h][k][i] = compute_infec_prof_inst(h, k, t_infec_i);
			LUT::cumul_profile[h][k][i] = compute_infec_prof_cumul(h, k, t_infec_i);
		}
		// NOTE : else is the case t_infec_i == t_infec_k, i. e. nothing to update
	}

	//=============================//
	//=====Susceptible Members=====//
	//=============================//

	t_infec_k = househ[h].get_data("end_followup");
	for (int k : househ[h].susceptible_set)
	{
		LUT::cumul_profile[h][i][k] = compute_infec_prof_cumul(h, i, t_infec_k);
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/* Initialize parameters and augmented variables */
void Model::common_init_model()
{
	// FIXME : the relative order of init_model_for_inference()
	// and init_augmented_variables() must be checked now

	//===== Initialize Lookup Tables (infectivity profiles) =====//

	// initialize the infectivity profile function
	LUT::infec_prof_func = Infec_Prof_Skew_Logistic(
		{{"tost_alpha", get_address("tost_alpha")},
		 {"tost_mu", get_address("tost_mu")},
		 {"tost_sigma", get_address("tost_sigma")},
		 {"tost_tau", get_address("tost_tau")}});
		 
	// clear the Lookup Tables (in case of multiple initializations)
	LUT::inst_profile.clear();
	LUT::cumul_profile.clear();

	//===== Initialize household features =====//

	// produce the empty matrix for the infectivity profiles look-up-tables
	for (int h = 0; h < num_househ; h++)
	{
		int h_size = househ[h].house_size;

		// produce the empty matrix for the infectivity profiles look-up-tables
		LUT::inst_profile.emplace_back(h_size, std::vector<double>(h_size, 0.));
		LUT::cumul_profile.emplace_back(h_size, std::vector<double>(h_size, 0.));
	}

	// find maximum study period within the household
	for (int h = 0; h < num_househ; h++)
	{
		double end_followup = 0;
		for (int i : househ[h].infectious_set)
		{
			end_followup = std::max(end_followup, househ[h].indiv[i].get_data("study_period"));
		}
		househ[h].assign_data("end_followup", end_followup);
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : hazards of infection  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Sum over all hazards of infection. */
double Model::hazards_of_infec(int h, int i, std::vector<squarematrix> &profiles)
{
	double output = 0.0;

	// Risk of infection from other household members (k)
	for (int k : househ[h].infectious_set)
	{
		if (k != i)
		{
			output += househ[h].get_contact_rate(k, i) * househ[h].indiv[k].get_group_rate("Infectivity") * profiles[h][k][i];
		}
	}

	return output * get_value("beta") * pow(get_value("median_size") / househ[h].house_size, get_value("delta"));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  /////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_house_loglike(int h)
{
	double loglike_house = 0.0;
	double t0 = get_infec_time(h, househ[h].first_case);
	double tf = househ[h].get_data("end_followup");
	double comm_risk_inst, comm_risk_cumul;

	comm_risk_cumul = get_value("alpha") * (tf - t0);

	// Contribution of healthy members
	for (int i : househ[h].susceptible_set)
	{
		loglike_house -= househ[h].get_group_rate(i, "Susceptibility") * (comm_risk_cumul + hazards_of_infec(h, i, LUT::cumul_profile));
	}

	// Contribution of symptomatic members
	for (int i : househ[h].infectious_set)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : augm_with_model_set)
		{
			loglike_house += househ[h].indiv[i].augm_model_log_pdf(augmID);
		}

		// Contribution of secondary cases
		if (i != househ[h].first_case)
		{
			comm_risk_inst = get_value("alpha");

			loglike_house += log(househ[h].get_group_rate(i, "Susceptibility") * (comm_risk_inst + hazards_of_infec(h, i, LUT::inst_profile)));
			loglike_house -= househ[h].get_group_rate(i, "Susceptibility") * (comm_risk_cumul + hazards_of_infec(h, i, LUT::cumul_profile));
		}
	}

	return loglike_house;
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Model : Methods (simulation)  ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate outbreak in a single household. */
void Model::simulate_single_house(int h)
{
	// note: incub periods are initialized according to option in augm file
	// note: also simulation are conditioned on the the first cases

	const double dt = 0.01;
	// FIXME : what if there is no first case?
	const double t0 = get_infec_time(h, househ[h].first_case);
	const double tf = househ[h].get_data("end_followup");

	//=================================================//
	//===== GENERATION 0 : by single introduction =====//
	//=================================================//

	// set all secondary cases to healthy
	for (int i = 0; i < househ[h].house_size; i++)
	{
		if (i != househ[h].first_case)
		{
			househ[h].change_status(i, SUSCEPTIBLE);
		}
	}

	//=============================================================//
	//===== GENERATION >=1 : infection by previous generation =====//
	//=============================================================//

	double esc_rate, comm_risk_inst;
	for (double t = dt + t0; (t < tf) && (househ[h].num_susceptible > 0); t += dt)
	{
		for (int i = 0; i < househ[h].house_size; i++)
		{
			if (househ[h].is_susceptible[i])
			{
				// set new status for computations
				househ[h].change_status(i, INFECTIOUS);
				househ[h].indiv[i].set_augm(
					"onset_time",
					t + uniform() * dt + househ[h].indiv[i].get_augm("incub_period"));
				augm_var_store_dependecies(h, i);

				comm_risk_inst = get_value("alpha");

				// refuse according to Bernoulli process
				esc_rate = househ[h].get_group_rate(i, "Susceptibility") * (comm_risk_inst + hazards_of_infec(h, i, LUT::inst_profile)) * dt;

				if (uniform() < exp(-esc_rate))
				{
					househ[h].change_status(i, SUSCEPTIBLE);
					augm_var_store_dependecies(h, i);
				}
			}
		}
	}
}