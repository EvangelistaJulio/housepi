// augm_debugger.cpp //

#include "augm_debugger.h"

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// Individual : Methods for Model ////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
double Model::get_infec_time(int i) const
{
	return indiv[i].get_augm("onset_time") - indiv[i].get_augm("incub_period");
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

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

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Model : Infectivity Profiles /////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Instanteneous infectivity profile. */
double Model::compute_infec_profile_inst(int i, double t_transm)
{
	double tost = t_transm - indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	double incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func.pdf(tost, incub_period);
}

/* Cumulative infectivity profile. */
double Model::compute_infec_profile_cumul(int i, double t_transm)
{
	double tost = t_transm - indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	double incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func.cdf(tost, incub_period);
}

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void Model::augm_var_store_dependecies(int i)
{
	int h = indiv[i].home;
	double t_infec_i = get_infec_time(i);
	double t_infec_k;

	//=====infectious members=====//
	for (int k : househ[h].infectious_list)
	{
		t_infec_k = get_infec_time(k);

		if (t_infec_i < t_infec_k)
		{ // `i` : infector, `k` : infectee
			inst_profile_LUT[from_to(i, k)] = compute_infec_profile_inst(i, t_infec_k);
			cumul_profile_LUT[from_to(i, k)] = compute_infec_profile_cumul(i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			inst_profile_LUT[from_to(k, i)] = compute_infec_profile_inst(k, t_infec_i);
			cumul_profile_LUT[from_to(k, i)] = compute_infec_profile_cumul(k, t_infec_i);
		}
		// NOTE : else is the case t_infec_i == t_infec_k, i. e. nothing to update
	}

	//=====susceptible members=====//
	t_infec_k = househ[h].get_data("end_followup");
	for (int k : househ[h].susceptibles_list)
	{
		cumul_profile_LUT[from_to(i, k)] = compute_infec_profile_cumul(i, t_infec_k);
	}
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/* Initialize parameters and augmented variables */
void Model::initialize_model()
{
	//-----initialize infectivity profiles-----//
	// dictionary with the address to the model parameters for the infectivity profile
	sdargs_ext infec_profile_args = {
		{"tost_alpha", get_address("tost_alpha")},
		{"tost_mu", get_address("tost_mu")},
		{"tost_sigma", get_address("tost_sigma")},
		{"tost_tau", get_address("tost_tau")}};

	// initialize the infectivity profile function of the model
	infec_profile_func = Infectivity_Profile_Function(infec_profile_args);

	// produce the empty matrix for the infectivity profiles look-up-tables
	inst_profile_LUT.assign(num_indiv * num_indiv, 0.0);
	cumul_profile_LUT.assign(num_indiv * num_indiv, 0.0);

	//-----find the first case-----//
	for (int h = 0; h < num_househ; h++)
	{
		// find maximum study period within the household
		double end_followup = 0;
		for (int i : househ[h].infectious_list)
		{
			end_followup = std::max(end_followup, indiv[i].get_data("study_period"));
		}
		househ[h].assign_data("end_followup", end_followup);

		// find the first case for each household
		update_first_case(h);
	}

	// initialize the dependencies of augmented variables
	// store relative hazard for each of the symptomatic and relative infectees
	for (int i = 0; i < num_indiv; i++)
	{
		augm_var_store_dependecies(i);
	}
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Sum over all hazards of infection. */
double Model::hazards_of_infec(std::vector<double> *profiles, int i)
{
	double output = 0.0;
	int h = indiv[i].home;

	// Risk of infection from other household members
	for (int k : househ[h].infectious_list)
	{
		if (k != i)
		{
			output += (*profiles)[from_to(k, i)];
		}
	}

	return output * get_value("beta") * pow(househ[h].house_size, -get_value("delta"));
}

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(int h)
{
	double loglike_house = 0.0;

	// Contribution of susceptible members
	for (int i : househ[h].susceptibles_list)
	{
		loglike_house -= hazards_of_infec(&cumul_profile_LUT, i);
	}

	// Contribution of symptomatic members
	for (int i : househ[h].infectious_list)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : augm_with_model_set)
		{
			loglike_house += indiv[i].augm_model_log_pdf(augmID);
		}

		// Contribution of secondary cases
		if (i != househ[h].first_case)
		{
			loglike_house += log(hazards_of_infec(&inst_profile_LUT, i));
			loglike_house -= hazards_of_infec(&cumul_profile_LUT, i);
		}
	}

	return loglike_house;
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

/* Poisson point process to simulate outbreak in a single household. */
void Model::simulate_single_house(int h)
{
	const double dt = 0.01;
	// note: incub period of all members are already initialized
	// note: first cases are already initialized
	const double t0 = get_infec_time(househ[h].first_case);
	const double tf = househ[h].get_data("end_followup");

	//-----set all secondary cases to healthy-----//
	std::vector<int> store_old_list = househ[h].infectious_list;
	for (int i : store_old_list)
	{
		if (i != househ[h].first_case)
		{
			househ[h].change_status(i, SUSCEPTIBLE);
		}
	}

	//=====simulate infection within household=====//
	double esc_proba;
	for (double t = dt + t0; (t < tf) && !househ[h].susceptibles_list.empty(); t += dt)
	{
		std::vector<int> store_old_list = househ[h].susceptibles_list;
		for (int i : store_old_list)
		{
			// set new status for computations
			househ[h].change_status(i, INFECTIOUS);
			indiv[i].set_augm("onset_time", t + uniform() * dt + indiv[i].get_augm("incub_period"));
			augm_var_store_dependecies(i);

			// refuse according to Bernoulli process
			esc_proba = exp(-hazards_of_infec(&inst_profile_LUT, i) * dt);

			if (uniform() < esc_proba)
			{
				househ[h].change_status(i, SUSCEPTIBLE);
				augm_var_store_dependecies(i);
			}
		}
	}
}

/* Output formatter */
void Model::format_augmented_data()
{
	for (int i = 0; i < num_indiv; i++)
	{
		if (indiv[i].get_status() == SUSCEPTIBLE)
		{
			indiv[i].set_augm("onset_time", 1000);
			indiv[i].set_augm("incub_period", 0);
		}
	}
}