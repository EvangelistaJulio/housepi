// user_model.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  * *** * AUGMENTED DEBUGGER * *** *  //////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "user_model.h"

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////// Infectivity_Profile_Function : probability distribution function /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Probability distribution function for the Infectivity Profile. */
double Infectivity_Profile_Function::pdf(double tost, double t_incub)
{
	double a = get_arg(0); // shape argument
	double m = get_arg(1); // location argument
	double s = get_arg(2); // shape argument
	double T = get_arg(3); // mean incubation period argument

	double norm_pdf = (s / a) * (1. + (t_incub / T - 1.) * pow(1 + exp(m / s), -a) - (t_incub / T) * (pow(1. + exp((m + T) / s), -a)));

	if (tost < 0)
	{
		tost *= (T / t_incub);
	}
	return exp((m - tost) / s) * pow(1. + exp((m - tost) / s), -1. - a) / norm_pdf;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// Infectivity_Profile_Function : cumulative distribution function /////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Cumulative distribution function forthe Infectivity Profile. */
double Infectivity_Profile_Function::cdf(double tost, double t_incub)
{
	double a = get_arg(0); // shape argument
	double m = get_arg(1); // location argument
	double s = get_arg(2); // shape argument
	double T = get_arg(3); // mean incubation period argument

	double norm_cdf = (1. + (t_incub / T - 1.) * pow(1 + exp(m / s), -a) - (t_incub / T) * (pow(1. + exp((m + T) / s), -a)));

	if (tost < 0)
	{
		// Integrate from -t_incub to tost
		return (t_incub / T) * (pow(1. + exp((m - tost * T / t_incub) / s), -a) - pow(1. + exp((m + T) / s), -a)) / norm_cdf;
	}
	else
	{
		// Integrate from -t_incub to 0
		double tmp = (t_incub / T) * (pow(1. + exp(m / s), -a) - pow(1. + exp((m + T) / s), -a));
		// Integrate from 0 to tost
		return (tmp + pow(1. + exp((m - tost) / s), -a) - pow(1. + exp(m / s), -a)) / norm_cdf;
	}
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  Model : get infection time  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
double Model::get_infec_time(int h, int i) const
{
	return househ[h].indiv[i].get_augm("onset_time") - househ[h].indiv[i].get_augm("incub_period");
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Model : compute infectivity profile instantaneous ///////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Instanteneous infectivity profile. */
double Model::compute_infec_profile_inst(int h, int i, double t_transm)
{
	// tost : time from onset of symptoms to transmission
	double tost = t_transm - househ[h].indiv[i].get_augm("onset_time");
	double incub_period = househ[h].indiv[i].get_augm("incub_period");
	return infec_profile_func.pdf(tost, incub_period);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Model : compute infectivity profile cumulative ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Cumulative infectivity profile. */
double Model::compute_infec_profile_cumul(int h, int i, double t_transm)
{
	double tost = t_transm - househ[h].indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	double incub_period = househ[h].indiv[i].get_augm("incub_period");
	return infec_profile_func.cdf(tost, incub_period);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// Model : augmented variables store dependencies //////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void Model::augm_var_store_dependecies(int h, int i)
{
	double t_infec_i = get_infec_time(h, i);
	double t_infec_k;

	//============================//
	//=====Infectious Members=====//
	//============================//

	for (int k : househ[h].infectious_list)
	{
		t_infec_k = get_infec_time(h, k);

		if (t_infec_i < t_infec_k)
		{ // `i` : infector, `k` : infectee
			inst_profile_LUT[h][i][k] = compute_infec_profile_inst(h, i, t_infec_k);
			cumul_profile_LUT[h][i][k] = compute_infec_profile_cumul(h, i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			inst_profile_LUT[h][k][i] = compute_infec_profile_inst(h, k, t_infec_i);
			cumul_profile_LUT[h][k][i] = compute_infec_profile_cumul(h, k, t_infec_i);
		}
		// NOTE : else is the case t_infec_i == t_infec_k, i. e. nothing to update
	}

	//=============================//
	//=====Susceptible Members=====//
	//=============================//

	t_infec_k = househ[h].get_data("end_followup");
	for (int k : househ[h].susceptibles_list)
	{
		cumul_profile_LUT[h][i][k] = compute_infec_profile_cumul(h, i, t_infec_k);
	}
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  Model : initialize model  ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize parameters and augmented variables */
void Model::initialize_model()
{
	//=====Initialize infectivity profiles=====//
	// dictionary with the address to the model parameters for the infectivity profile
	orderargs_ptr infec_profile_args = {
		{"tost_alpha", get_param_pointer("tost_alpha")},
		{"tost_mu", get_param_pointer("tost_mu")},
		{"tost_sigma", get_param_pointer("tost_sigma")},
		{"tost_tau", get_param_pointer("tost_tau")}};

	// initialize the infectivity profile function of the model
	infec_profile_func = Infectivity_Profile_Function(infec_profile_args);

	//=====Initialize household features=====//
	for (int h = 0; h < num_househ; h++)
	{	
		// produce the empty matrix for the infectivity profiles look-up-tables
		int h_size = househ[h].house_size;
		inst_profile_LUT.emplace_back(h_size, std::vector<double>(h_size, 0.0));
		cumul_profile_LUT.emplace_back(h_size, std::vector<double>(h_size, 0.0));

		// find maximum study period within the household
		double end_followup = 0;
		for (int i : househ[h].infectious_list)
		{
			end_followup = std::max(end_followup, househ[h].indiv[i].get_data("study_period"));
		}
		househ[h].assign_data("end_followup", end_followup);

		// find the first case for each household
		update_first_case(h);
	}

	//=====Initialize the dependencies of augmented variables=====//
	// store relative hazard for each of the symptomatic and relative infectees
	for (int h = 0; h < num_househ; h++)
	{
		for (int i = 0; i < househ[h].house_size; i++)
		{
			augm_var_store_dependecies(h, i);
		}
	}
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : hazards of infection  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Sum over all hazards of infection. */
double Model::force_of_infec(int h, int i, std::vector<squarematrix> &profiles)
{
	double output = 0.0;

	// Risk of infection from other household members
	for (int k : househ[h].infectious_list)
	{
		if (k != i)
		{
			output += profiles[h][k][i];
		}
	}

	return output * get_param("beta") * pow(househ[h].house_size, -get_param("delta"));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house log-likelihood  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(int h)
{
	double loglike_house = 0.0;

	//=============================//
	//=====Susceptible Members=====//
	//=============================//

	for (int i : househ[h].susceptibles_list)
	{
		loglike_house -= force_of_infec(h, i, cumul_profile_LUT);
	}

	//============================//
	//=====Infectious Members=====//
	//============================//

	for (int i : househ[h].infectious_list)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : augm_with_model_set)
		{
			loglike_house += househ[h].indiv[i].augm_model_log_pdf(augmID);
		}

		// Contribution of secondary cases
		if (i != househ[h].first_case)
		{
			loglike_house += log(force_of_infec(h, i, inst_profile_LUT));
			loglike_house -= force_of_infec(h, i, cumul_profile_LUT);
		}
	}

	return loglike_house;
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Model : simulate single house  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate outbreak in a single household. */
void Model::simulate_single_house(int h)
{
	// note: incub periods are initialized according to option in augm file
	// note: first cases are already initialized

	const double dt = 0.01;
	const double t0 = get_infec_time(h, househ[h].first_case);
	const double tf = househ[h].get_data("end_followup");

	//=====Set all secondary cases to healthy=====//
	std::vector<int> store_old_list = househ[h].infectious_list;
	for (int i : store_old_list)
	{
		if (i != househ[h].first_case)
		{
			househ[h].change_status(i, SUSCEPTIBLE);
		}
	}

	//=====Simulate infection within household=====//
	double esc_proba;
	for (double t = dt + t0; (t < tf) && !househ[h].susceptibles_list.empty(); t += dt)
	{
		std::vector<int> store_old_list = househ[h].susceptibles_list;
		for (int i : store_old_list)
		{
			// set new status for computations
			househ[h].change_status(i, INFECTIOUS);
			househ[h].indiv[i].set_augm("onset_time", t + uniform() * dt + househ[h].indiv[i].get_augm("incub_period"));
			augm_var_store_dependecies(h, i);

			// refuse according to Bernoulli process
			esc_proba = exp(-force_of_infec(h, i, inst_profile_LUT) * dt);

			if (uniform() < esc_proba)
			{
				househ[h].change_status(i, SUSCEPTIBLE);
				augm_var_store_dependecies(h, i);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////  Model : format augmented data  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Format the augmented variables before printing them to output. */
void Model::format_augmented_data()
{
	for (int h = 0; h < num_househ; h++)
	{
		for (int i = 0; i < househ[h].house_size; i++)
		{
			if (househ[h].indiv[i].get_status() == SUSCEPTIBLE)
			{
				househ[h].indiv[i].set_augm("onset_time", 1000);
				househ[h].indiv[i].set_augm("incub_period", 0);
			}
		}
	}
}

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//