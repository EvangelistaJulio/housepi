// model-simple.cpp //

#include "model.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////  Extra Files ////////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Read Extra File */
void build_extra(string data_pathtofile)
{
	UNUSED(data_pathtofile);
}

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
// compute likelihood straight
void Model::store_augm_var_dependecies(size_t i)
{
	UNUSED(i);
}

/* Assign the contact patterns between a couple of individuals. */
//
void Model::assign_rate_per_contact(size_t i, size_t k)
{
	UNUSED(i);
	UNUSED(k);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////// Model : Augmented Variables //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* List of Model Quantities */
// FIXME (temporary) : this should be automatic
void Model::fill_model_lists()
{
	// model group rates
	group_rate_list.push_back("Infectivity");
	group_rate_list.push_back("Susceptibility");
}

/* Initialize augmented variables (MANDATORY) */
// FIXME: this needs to become more user friendly
void Model::assign_augm_to_individual_(Individual *this_indiv, std::mt19937_64 &rng)
{
	int status = this_indiv->get_status();

	//=====incub_period=====//
	{
		std::string augmID = "incub_period";
		double init_value = 0.;

		std::string model_name = (status == HEALTHY) ? UNDEFINED : "gamma";
		kwargs_ext model_args = {
			{"shape", get_address("incub_shape")},
			{"scale", get_address("incub_scale")}};

		std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "normal";
		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", get_value("incub_mcmc_sd")}};
		double metropolis_zero = INFINITY; // from old value

		/*
		// FIXME : prior requires a checkup
		std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "prior";
		kwargs metropolis_args = {};
		*/

		kwargs range = {
			{"low", 0.0},
			{"up", get_value("incub_max_value")}};

		Augmented augm_var(
			init_value, model_name, model_args,
			metropolis_strategy, metropolis_args, metropolis_zero, range);
		(status == HEALTHY) ? this_indiv->add_augm(augmID, augm_var) : this_indiv->add_augm(augmID, augm_var, rng);

		// automatic insertion to sets
		augm_set.insert(augmID);
		if (model_name != UNDEFINED)
		{
			augm_with_model_set.insert(augmID);
		}
		if (metropolis_strategy != UNDEFINED)
		{
			augm_to_infer_set.insert(augmID);
		}
	}

	//=====onset_time=====//
	{
		std::string augmID = "onset_time";
		double init_value = (status == HEALTHY) ? this_indiv->get_data("study_period") : this_indiv->get_data("onset_time");

		std::string model_name = UNDEFINED;
		kwargs_ext model_args = {};

		std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "normal";

		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", get_value("onset_err")}};

		double metropolis_zero = (status == HEALTHY) ? INFINITY : init_value;
		kwargs range = {
			{"low", get_value("onset_low")},
			{"up", this_indiv->get_data("study_period")}}; // FIXME onset_time == study_period!!!

		Augmented augm_var(
			init_value, model_name, model_args,
			metropolis_strategy, metropolis_args, metropolis_zero, range);
		this_indiv->add_augm(augmID, augm_var);

		// automatic insertion to sets
		augm_set.insert(augmID);
		if (model_name != UNDEFINED)
		{
			augm_with_model_set.insert(augmID);
		}
		if (metropolis_strategy != UNDEFINED)
		{
			augm_to_infer_set.insert(augmID);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// Individual : Methods for Model ////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
timetype Individual::get_infec_time() const
{
	return get_augm("onset_time") - get_augm("incub_period");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// Infectivity_Profile_Function //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Choice of the Infectivity Profile Function.
The infectivity profile is considered equal to the Time from Onset of Symptoms to Transmission `tost`
>	tost = <time infection infectee> - <time onset infector>
Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
where the pre-symptomaticvalues are rescaled as a function of the incubation time `t_incub`.
The arguments to provide are :
	0) "tost_alpha" : shape;
	1) "tost_mu" : location;
	2) "tost_sigma" : scale;
	3) "tost_tau" : mean incubation time;

Side note: symtoms can occur after infection (pre-symptomatic transmission), in other words `tost` can be negative.

----
Ref:
----
Ferretti, L. et al. The timing of COVID-19 transmission. 2020.09.04.20188516
Preprint at https://doi.org/10.1101/2020.09.04.20188516 (2020).*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Model : Infectivity Profiles /////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the infectivity profiles. */
void Model::assign_infec_profile_function()
{
	// dictionary with the address to the model parameters for the infectivity profile
	kwargs_ext infec_profile_args = {
		{"tost_alpha", get_address("tost_alpha")},
		{"tost_mu", get_address("tost_mu")},
		{"tost_sigma", get_address("tost_sigma")},
		{"tost_tau", get_address("tost_tau")}};

	// initialize the infectivity profile function of the model
	infec_profile_func = new Infectivity_Profile_Function(infec_profile_args);
}

/* Instanteneous infectivity profile. */
double Model::compute_infec_profile_inst(size_t i, timetype t_transm)
{
	timetype tost = t_transm - indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	timetype incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func->pdf(tost, incub_period);
}

/* Cumulative infectivity profile. */
double Model::compute_infec_profile_cumul(size_t i, timetype t_transm)
{
	timetype tost = t_transm - indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	timetype incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func->cdf(tost, incub_period);
}

/* Probability distribution function for the Infectivity Profile. */
double Infectivity_Profile_Function::pdf(timetype tost, timetype t_incub)
{
	double a = get_arg(0); // shape argument
	double m = get_arg(1); // location argument
	double s = get_arg(2); // shape argument
	double T = get_arg(3); // mean incubation period argument

	double norm_pdf = (s / a) * (1. + (t_incub / T - 1.) * std::pow(1 + std::exp(m / s), -a) - (t_incub / T) * (std::pow(1. + std::exp((m + T) / s), -a)));

	if (tost < 0)
	{
		tost *= (T / t_incub);
	}
	return std::exp((m - tost) / s) * std::pow(1. + std::exp((m - tost) / s), -1. - a) / norm_pdf;
}

/* Cumulative distribution function forthe Infectivity Profile. */
double Infectivity_Profile_Function::cdf(timetype tost, timetype t_incub)
{
	double a = get_arg(0); // shape argument
	double m = get_arg(1); // location argument
	double s = get_arg(2); // shape argument
	double T = get_arg(3); // mean incubation period argument

	double norm_cdf = (1. + (t_incub / T - 1.) * std::pow(1 + std::exp(m / s), -a) - (t_incub / T) * (std::pow(1. + std::exp((m + T) / s), -a)));

	if (tost < 0)
	{
		// Integrate from -t_incub to tost
		return (t_incub / T) * (std::pow(1. + std::exp((m - tost * T / t_incub) / s), -a) - std::pow(1. + std::exp((m + T) / s), -a)) / norm_cdf;
	}
	else
	{
		// Integrate from -t_incub to 0
		double tmp = (t_incub / T) * (std::pow(1. + std::exp(m / s), -a) - std::pow(1. + std::exp((m + T) / s), -a));
		// Integrate from 0 to tost
		return (tmp + std::pow(1. + std::exp((m - tost) / s), -a) - std::pow(1. + std::exp(m / s), -a)) / norm_cdf;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////// Model : associated rates  ////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the parameters controlling infectivity for each individual */
// assign to
void Model::assign_rate_per_group_to_individual(Individual *this_indiv)
{
	//=====rate 0 : Infectivity=====//
	if ((this_indiv->get_data("recov") == false) && (this_indiv->get_data("vacc_status") == 0))
	{ // unrecovered unvaccinated
		this_indiv->assign_rate_per_group(0, "i_naive");
	}
	else
	{ // vaccinated o recovered : Infectivity [REFERENCE]
		this_indiv->assign_rate_per_group(0, REFERENCE);
	}

	//=====rate 1 : Susceptiblity=====//
	if (this_indiv->get_data("vacc_status") == 0)
	{ // naive : Susceptiblity [REFERENCE]
		this_indiv->assign_rate_per_group(1, REFERENCE);
	}
	else
	{ // 1 dose
		this_indiv->assign_rate_per_group(1, "s_vacc");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_loglike(size_t h)
{
	double hazard_inst = 0.0, hazard_cumul = 0.0;
	double loglike_house = 0.0;
	size_t housesize = househ[h].size;
	timetype t0 = househ[h].start_followup;
	timetype t_infec_i, t_infec_k;

	//=====loop over all healthy residents=====//
	t_infec_i = househ[h].end_followup;
	for (size_t i : househ[h].susceptibles_list)
	{
		// loop over all possible infectors
		hazard_cumul = 0.0;
		for (size_t k : househ[h].infectors_list)
		{
			hazard_cumul += indiv[k].model_rate(0) * compute_infec_profile_cumul(k, t_infec_i);
		}
		loglike_house -= indiv[i].model_rate(1) * (get_value("alpha") * (t_infec_i - t0) + get_value("beta") * std::pow(get_value("median_size") / housesize, get_value("delta")) * hazard_cumul);
	}

	//=====loop over all infectprs residents=====//
	for (size_t i : househ[h].infectors_list)
	{
		t_infec_i = indiv[i].get_infec_time();

		// Contribution of model associated to augmented variables
		for (std::string augmID : augm_with_model_set)
		{
			loglike_house += indiv[i].augm_model_log_pdf(augmID);
		}
		
		// Contribution of secondary cases
		if (i != househ[h].first_case)
		{
			// loop over all possible infectors
			hazard_inst = 0.0;
			hazard_cumul = 0.0;
			for (size_t k : househ[h].infectors_list)
			{
				t_infec_k = indiv[k].get_infec_time();

				if (t_infec_i < t_infec_k)
				{ // `i` : infector, `k` : infectee
					hazard_inst += compute_infec_profile_inst(i, t_infec_k);
					hazard_cumul += compute_infec_profile_cumul(i, t_infec_k);
				}
				else if (t_infec_k < t_infec_i)
				{ // `k` : infector, `i` : infectee
					hazard_inst += compute_infec_profile_inst(k, t_infec_i);
					hazard_cumul += compute_infec_profile_cumul(k, t_infec_i);
				}
				// NOTE : else is the case t_infec_i == t_infec_k, i. e. nothing to update
			}
		}
		loglike_house += std::log(indiv[i].model_rate(1) * (get_value("alpha") + get_value("beta") * pow(househ[h].size, -get_value("delta")) * hazard_inst));
		loglike_house -= indiv[i].model_rate(1) * (get_value("alpha") * (t_infec_i - t0) + get_value("beta") * pow(househ[h].size, -get_value("delta")) * hazard_cumul);
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

/* Poisson point process to simulate data. */
void Model::simulate_infection(timetype dt, std::mt19937_64 &rng)
{
	double proba, hazard_inst;
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	timetype end_followup;

	// simulate the infection for each household
	for (size_t h = 0; h < num_househ; h++)
	{
		// fix all incub period at the beginning (all others should already be initialized)
		for (size_t i : househ[h].susceptibles_list)
		{
			indiv[i].init_augm("incub_period", rng);
		}

		// set all secondary cases to healthy
		std::vector<size_t> copy_list = househ[h].infectors_list;
		for (size_t i : copy_list)
		{
			if (i != househ[h].first_case)
			{
				indiv[i].set_status(HEALTHY);
				househ[h].change_status(i, HEALTHY);
			}
		}

		//=====simulate infection within household=====//
		end_followup = househ[h].end_followup;
		for (timetype t = dt + househ[h].start_followup; ((t < end_followup) && !househ[h].susceptibles_list.empty()); t += dt)
		{
			std::vector<size_t> copy_list = househ[h].susceptibles_list;
			for (size_t i : copy_list)
			{
				househ[h].change_status(i, SYMPTOMATIC);
				indiv[i].set_status(SYMPTOMATIC);

				indiv[i].set_augm("onset_time", t + unif(rng) * dt + indiv[i].get_augm("incub_period"));

				// Contribution of secondary cases

				// loop over all possible infectors
				hazard_inst = 0.0;
				for (size_t k : househ[h].infectors_list)
				{

					hazard_inst += compute_infec_profile_inst(k, t);
				}
				// proba of infection
				proba = 1 - exp(-indiv[i].model_rate(1) * (get_value("alpha") + get_value("beta") * pow(househ[h].size, -get_value("delta")) * hazard_inst) * dt);

				// refuse
				if (unif(rng) >= proba)
				{
					househ[h].change_status(i, HEALTHY);
					indiv[i].set_status(HEALTHY);
				}
			}
		}
	}

	// Output formatter
	for (size_t i = 0; i < num_indiv; i++)
	{
		if (indiv[i].get_status() == SYMPTOMATIC)
		{
			if (indiv[i].get_augm("onset_time") > indiv[i].get_data("study_period"))
			{
				indiv[i].set_status(HEALTHY);
			}
			else
			{
				indiv[i].set_augm("onset_time", round(indiv[i].get_augm("onset_time")));
			}
		}

		if (indiv[i].get_status() == HEALTHY)
		{
			indiv[i].set_augm("onset_time", 1000);
			indiv[i].set_augm("incub_period", 0);
		}
	}
}