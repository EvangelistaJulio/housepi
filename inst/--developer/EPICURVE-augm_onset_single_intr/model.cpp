// model.cpp //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

#include "model.h"

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////  Extra Files ////////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

// use the EPI curve for the community contribution
std::vector<double> Date_EPI_Curve = {};
std::vector<double> Inst_EPI_Curve = {};
std::vector<double> Cumul_EPI_Curve = {};

/* Read EPI curve from csv file. */
void epi_curve_read_csv(std::string pathtofile)
{
	std::ifstream input(pathtofile.c_str());
	if (input.is_open())
	{
		// load instantaneous value of the epi curve
		std::string line;
		while (getline(input, line))
		{
			std::stringstream ss(line);
			std::string item;

			// Store information in variables
			if (std::getline(ss, item, ','))
				Date_EPI_Curve.push_back(std::stod(item));
			if (std::getline(ss, item, ','))
				Inst_EPI_Curve.push_back(std::stod(item));
		}

		// store the cumulative epi curve
		Cumul_EPI_Curve.push_back(Inst_EPI_Curve.front());
		for (size_t i = 1; i < Inst_EPI_Curve.size(); ++i)
			Cumul_EPI_Curve.push_back(Inst_EPI_Curve[i] + Cumul_EPI_Curve[i - 1]);
	}
	else
		tools::error("Impossible to read file: " + pathtofile + ".");
}

/* Linear interpolation by date. */
double interpolate_by_date(double x0, const std::vector<double> &y)
{
	const std::vector<double> &x = Date_EPI_Curve; // alias

	if (x0 <= x.front()) // saturate to lower value
		return y.front();
	else if (x0 >= x.back()) // saturate to upper value
		return y.back();

	// else calculate the interpolated y value
	// bisection
	size_t i = 0;			 // lower bound
	size_t u = x.size() - 1; // upper bound
	while (u - i > 1)
	{
		size_t mid = i + (u - i) / 2;
		if (x[mid] <= x0)
			i = mid;
		else
			u = mid;
	}
	return y[i] + (x0 - x[i]) * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
}

/* Interpolated value of the instantaneous epi curve. */
double epi_curve_inst(double t) { return interpolate_by_date(t, Inst_EPI_Curve); }

/* Interpolated value of the cumulative epi curve. */
double epi_curve_cumul(double t) { return interpolate_by_date(t, Cumul_EPI_Curve); }

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
	return house_vector[h].members[i].get_augm("onset_time") - house_vector[h].members[i].get_augm("incub_period");
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
	double tost = t_transm - house_vector[h].members[i].get_augm("onset_time");
	double incub_period = house_vector[h].members[i].get_augm("incub_period");
	return infec_profile_func.pdf(tost, incub_period);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////// Model : compute infectivity profile cumulative ////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Cumulative infectivity profile. */
double Model::compute_infec_profile_cumul(int h, int i, double t_transm)
{
	double tost = t_transm - house_vector[h].members[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	double incub_period = house_vector[h].members[i].get_augm("incub_period");
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

	for (int k : house_vector[h].get_infected_status_hhindexes)
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

	t_infec_k = house_vector[h].data("end_followup");
	for (int k : house_vector[h].susceptibles_list)
	{
		cumul_profile_LUT[h][i][k] = compute_infec_profile_cumul(h, i, t_infec_k);
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
	//=====Initialize infectivity profiles=====//
	// dictionary with the address to the model parameters for the infectivity profile
	orderargs_ptr infec_profile_args = {
		{"tost_alpha", param_pointer("tost_alpha")},
		{"tost_mu", param_pointer("tost_mu")},
		{"tost_sigma", param_pointer("tost_sigma")},
		{"tost_tau", param_pointer("tost_tau")}};

	// initialize the infectivity profile function of the model
	infec_profile_func = Infectivity_Profile_Function(infec_profile_args);

	//=====Initialize household features=====//
	for (int h = 0; h < num_house; h++)
	{
		// produce the empty matrix for the infectivity profiles look-up-tables
		int h_size = house_vector[h].house_size;
		inst_profile_LUT.emplace_back(h_size, std::vector<double>(h_size, 0.0));
		cumul_profile_LUT.emplace_back(h_size, std::vector<double>(h_size, 0.0));

		// find maximum study period within the household
		double end_followup = 0;
		for (int i : house_vector[h].get_infected_status_hhindexes)
		{
			end_followup = std::max(end_followup, house_vector[h].members[i].data("study_period"));
		}
		house_vector[h].assign_data("end_followup", end_followup);

		// find the first case for each household
		update_first_case(h);
	}

	//=====Initialize the dependencies of augmented variables=====//
	// store relative hazard for each of the symptomatic and relative infectees
	for (int h = 0; h < num_house; h++)
	{
		for (int i = 0; i < house_vector[h].house_size; i++)
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
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

void Model::assign_group_rate_to_individual(Individual &this_indiv)
{
	//=====rate 0 : Infectivity=====//
	if ((this_indiv.data("recov") == false) && (this_indiv.data("vacc_status") == 0))
	{ // unrecovered unvaccinated
		this_indiv.assign_group_rate("Infectivity", "i_naive");
	}
	else
	{ // vaccinated o recovered : Infectivity [REFERENCE]
		this_indiv.assign_group_rate("Infectivity", REFERENCE);
	}

	//=====rate 1 : Susceptibility=====//
	if (this_indiv.data("vacc_status") == 0)
	{ // naive : Susceptibility [REFERENCE]
		this_indiv.assign_group_rate("Susceptibility", REFERENCE);
	}
	else
	{ // 1 dose
		this_indiv.assign_group_rate("Susceptibility", "s_vacc");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  Model : assign rate_per group to individual  //////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_contact_rate_to_household(int h, int i, int k)
{
	// from `i` -> to `k`
	// children if younger than 12 y.o., adults otherwise
	bool i_is_child = (house_vector[h].members[i].data("age") < param("child_age"));
	bool i_is_adult = (i_is_child == false);

	bool k_is_child = (house_vector[h].members[k].data("age") < param("child_age"));
	bool k_is_adult = (k_is_child == false);

	if (i_is_child && k_is_child)
	{ // Child to Child
		house_vector[h].assign_contact_rate(i, k, "c_kid_to_kid");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		house_vector[h].assign_contact_rate(i, k, "c_kid_to_adu");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		house_vector[h].assign_contact_rate(i, k, "c_adu_to_kid");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		house_vector[h].assign_contact_rate(i, k, REFERENCE);
	}
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Model : hazards of infection  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Sum over all hazards of infection. */
double Model::force_of_infec(int h, int i, std::vector<matrix> &profiles)
{
	double output = 0.0;

	// Risk of infection from other household members (k)
	for (int k : house_vector[h].get_infected_status_hhindexes)
	{
		if (k != i)
		{
			output += house_vector[h].contact_rate(k, i) * house_vector[h].members[k].get_group_rate("Infectivity") * profiles[h][k][i];
		}
	}

	return output * param("beta") * pow(param("median_size") / house_vector[h].house_size, param("delta"));
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  Model : compute single house_vector log-likelihood  ///////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(int h)
{
	double loglike_house = 0.0;
	double t0 = get_infec_time(h, house_vector[h].first_case);
	double tf = house_vector[h].data("end_followup");
	double comm_risk_inst, comm_risk_cumul;

	// #ifdef TEST_EPI_CURVE
	comm_risk_cumul = get_value("alpha") * (interpolate_by_idx(tf, Cumul_EPI_Curve) - interpolate_by_idx(t0, Cumul_EPI_Curve));

	// Contribution of healthy members
	for (int i : house_vector[h].susceptibles_list)
	{
		loglike_house -= house_vector[h].members[i].get_group_rate("Susceptibility") * (comm_risk_cumul + force_of_infec(h, i, cumul_profile_LUT));
	}

	// Contribution of symptomatic members
	for (int i : house_vector[h].get_infected_status_hhindexes)
	{
		// Contribution of model associated to augmented variables
		for (std::string augm_key : augm_with_model_set)
		{
			loglike_house += house_vector[h].members[i].augm_model_logpdf(augm_key);
		}

		// Contribution of secondary cases
		if (i != house_vector[h].first_case)
		{
			// #ifdef TEST_EPI_CURVE
			comm_risk_inst = get_value("alpha") * interpolate_by_idx(get_infec_time(i), Inst_EPI_Curve);

			loglike_house += log(house_vector[h].members[i].get_group_rate("Susceptibility") * (comm_risk_inst + force_of_infec(h, i, inst_profile_LUT)));
			loglike_house -= house_vector[h].members[i].get_group_rate("Susceptibility") * (comm_risk_cumul + force_of_infec(h, i, cumul_profile_LUT));
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
void Model::simulate_single_household(int h)
{
	// note: incub periods are initialized according to option in augm file
	// note: first cases are already initialized

	const double dt = 0.01;
	const double t0 = get_infec_time(h, house_vector[h].first_case);
	const double tf = house_vector[h].data("end_followup");

	//-----set all secondary cases to healthy-----//
	std::vector<int> store_old_list = house_vector[h].get_infected_status_hhindexes;
	for (int i : store_old_list)
	{
		if (i != house_vector[h].first_case)
		{
			house_vector[h].change_status(i, HEALTHY_STATUS);
		}
	}

	//=====simulate infection within household=====//
	double esc_proba, comm_risk_inst;
	for (double t = dt + t0; (t < tf) && !house_vector[h].susceptibles_list.empty(); t += dt)
	{
		std::vector<int> store_old_list = house_vector[h].susceptibles_list;
		for (int i : store_old_list)
		{
			// set new status for computations
			house_vector[h].change_status(i, INFECTED_STATUS);
			house_vector[h].members[i].set_augm(
				"onset_time",
				t + uniform() * dt + house_vector[h].members[i].get_augm("incub_period"));
			augm_var_store_dependecies(h, i);

			// #ifdef TEST_EPI_CURVE
			comm_risk_inst = get_value("alpha") * interpolate_by_idx(get_infec_time(i), Inst_EPI_Curve);

			// refuse according to Bernoulli process
			esc_proba = exp(-house_vector[h].members[i].get_group_rate("Susceptibility") * (comm_risk_inst + force_of_infec(h, i, inst_profile_LUT)) * dt);

			if (uniform() < esc_proba)
			{
				house_vector[h].change_status(i, HEALTHY_STATUS);
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
	for (int h = 0; h < num_house; h++)
	{
		for (int i = 0; i < house_vector[h].house_size; i++)
		{
			if (house_vector[h].members[i].get_status() == HEALTHY_STATUS)
			{
				house_vector[h].members[i].set_augm("onset_time", 1000);
				house_vector[h].members[i].set_augm("incub_period", 0);
			}
		}
	}
}