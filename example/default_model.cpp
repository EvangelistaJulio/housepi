// default_model.cpp //

#include "default_model.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// Individual : Methods for Model ////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
double Individual::get_infec_time() const
{
	return get_augm("onset_time") - get_augm("incub_period");
}

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// Infectivity_Profile_Function //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Define the Infectivity_Profile_Function class */
class Infectivity_Profile_Function : public Base_Function
{
public:
	// Constructors
	Infectivity_Profile_Function(sdargs new_args) : Base_Function(
		new_args,
		{"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}
		) {}
	Infectivity_Profile_Function(sdargs_ext new_args_ext) : Base_Function(
		new_args_ext,
		{"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}
		) {}
	~Infectivity_Profile_Function(){};

	double pdf(double tost, double t_incub);
	double cdf(double tost, double t_incub);
};

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

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Model : Infectivity Profiles /////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the infectivity profiles. */
void Model::assign_infec_profile_function()
{
	// dictionary with the address to the model parameters for the infectivity profile
	sdargs_ext infec_profile_args = {
		{"tost_alpha", get_address("tost_alpha")},
		{"tost_mu", get_address("tost_mu")},
		{"tost_sigma", get_address("tost_sigma")},
		{"tost_tau", get_address("tost_tau")}};

	// initialize the infectivity profile function of the model
	infec_profile_func = new Infectivity_Profile_Function(infec_profile_args);
}

/* Instanteneous infectivity profile. */
double Model::compute_infec_profile_inst(size_t i, double t_transm)
{
	double tost = t_transm - indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	double incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func->pdf(tost, incub_period);
}

/* Cumulative infectivity profile. */
double Model::compute_infec_profile_cumul(size_t i, double t_transm)
{
	double tost = t_transm - indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	double incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func->cdf(tost, incub_period);
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
    assign_infec_profile_function(); // initialize infectivity profiles

    // find the first case for each household
    for (size_t h = 0; h < num_househ; h++)
    {
        househ[h].assign_end_followup();
        househ[h].update_first_case();
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
		this_indiv->assign_rate_per_group("Infectivity", "i_naive");
	}
	else
	{ // vaccinated o recovered : Infectivity [REFERENCE]
		this_indiv->assign_rate_per_group("Infectivity", REFERENCE);
	}

	//=====rate 1 : Susceptiblity=====//
	if (this_indiv->get_data("vacc_status") == 0)
	{ // naive : Susceptiblity [REFERENCE]
		this_indiv->assign_rate_per_group("Susceptibility", REFERENCE);
	}
	else
	{ // 1 dose
		this_indiv->assign_rate_per_group("Susceptibility", "s_vacc");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(size_t h)
{
	double hazard_inst = 0.0, hazard_cumul = 0.0;
	double loglike_house = 0.0;
	size_t housesize = househ[h].size;
	double t0 = househ[h].start_followup;
	double t_infec_i, t_infec_k;

	//=====loop over all healthy residents=====//
	t_infec_i = househ[h].end_followup;
	for (size_t i : househ[h].susceptibles_list)
	{
		// loop over all possible infectors
		hazard_cumul = 0.0;
		for (size_t k : househ[h].infectors_list)
		{
			hazard_cumul += indiv[k].model_rate("Infectivity") * compute_infec_profile_cumul(k, t_infec_i);
		}
		loglike_house -= indiv[i].model_rate("Susceptibility") * (get_value("alpha") * (t_infec_i - t0) + get_value("beta") * std::pow(get_value("median_size") / housesize, get_value("delta")) * hazard_cumul);
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
		loglike_house += std::log(indiv[i].model_rate("Susceptibility") * (get_value("alpha") + get_value("beta") * pow(househ[h].size, -get_value("delta")) * hazard_inst));
		loglike_house -= indiv[i].model_rate("Susceptibility") * (get_value("alpha") * (t_infec_i - t0) + get_value("beta") * pow(househ[h].size, -get_value("delta")) * hazard_cumul);
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
void Model::simulate_infection(double dt)
{
	double proba, hazard_inst;
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	double end_followup;

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
		for (double t = dt + househ[h].start_followup; ((t < end_followup) && !househ[h].susceptibles_list.empty()); t += dt)
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
				proba = 1 - exp(-indiv[i].model_rate("Susceptibility") * (get_value("alpha") + get_value("beta") * pow(househ[h].size, -get_value("delta")) * hazard_inst) * dt);

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