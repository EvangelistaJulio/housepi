// halloween.cpp //

#include "halloween.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////////  Extra Files ////////////////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

// use the EPI curve for the community contribution
std::vector<double> Cumul_EPI_Curve;
std::vector<double> Inst_EPI_Curve;

/* Read Extra File */
void build_extra(std::string data_pathtofile)
{
	double date, inst_val, cumul_val;

	std::ifstream input(data_pathtofile.c_str());
	if (input.is_open())
	{
		if (Verbose)
		{
			std::cout << "\n================EXTRA=================" << std::endl;
			std::cout << "Data file: " << data_pathtofile << std::endl;
		}

		// Create a std::stringstream of the current line
		for (std::string line; getline(input, line);)
		{
			std::istringstream iss(line);

			// Store information in variables
			iss >> date >> inst_val >> cumul_val;

			Inst_EPI_Curve.push_back(inst_val);
			Cumul_EPI_Curve.push_back(cumul_val);
		}
	}
	else
	{
		std::string err_text = "Impossible to read file: " + data_pathtofile + ".";
		exit_with_failure_error(err_text);
	}
}

/* Linear interpolation for an indexed std::vector. */
double lin_idx_interp(double x_value, const std::vector<double> &y)
{
	std::size_t x_0 = (int)floor(x_value);
	double output = 0.0;

	if (x_0 <= 0)
	{ // Saturate to lower value
		output = y[0];
	}
	else if (x_0 >= y.size() - 1)
	{ // Saturate to upper value
		output = y.back();
	}
	else
	{ // Calculate the interpolated y value
		output = y[x_0] + (x_value - x_0) * (y[x_0 + 1] - y[x_0]);
	}

	return output;
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

/* Define the Infectivity_Profile_Function class */
class Infectivity_Profile_Function : public Base_Function
{
public:
	// Constructors
	Infectivity_Profile_Function(sdargs new_args) : Base_Function(new_args, {"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}) {}
	Infectivity_Profile_Function(sdargs_ext new_args_ext) : Base_Function(new_args_ext, {"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}) {}
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

// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//
// -----------------------------------------------------------------------------------------------------------//

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

	// produce the empty matrix for the infectivity profiles look-up-tables
	inst_profile_lookuptab = new double[num_indiv * num_indiv];
	cumul_profile_lookuptab = new double[num_indiv * num_indiv];
}

/* Instanteneous infectivity profile. */
double Model::compute_infec_profile_inst(std::size_t i, double t_transm)
{
	double tost = t_transm - indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	double incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func->pdf(tost, incub_period);
}

/* Cumulative infectivity profile. */
double Model::compute_infec_profile_cumul(std::size_t i, double t_transm)
{
	double tost = t_transm - indiv[i].get_augm("onset_time"); // time from onset of symptoms to transmission
	double incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func->cdf(tost, incub_period);
}

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void Model::store_augm_var_dependecies(std::size_t i)
{
	std::size_t h = indiv[i].home;
	double t_infec_i = get_infec_time(i);
	double t_infec_k;

	//=====loop over all symptomatic residents=====//
	for (std::size_t k : househ[h].infectious_list)
	{
		t_infec_k = get_infec_time(k);

		if (t_infec_i < t_infec_k)
		{ // `i` : infector, `k` : infectee
			inst_profile_lookuptab[from_to(i, k)] = compute_infec_profile_inst(i, t_infec_k);
			cumul_profile_lookuptab[from_to(i, k)] = compute_infec_profile_cumul(i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			inst_profile_lookuptab[from_to(k, i)] = compute_infec_profile_inst(k, t_infec_i);
			cumul_profile_lookuptab[from_to(k, i)] = compute_infec_profile_cumul(k, t_infec_i);
		}
		// NOTE : else is the case t_infec_i == t_infec_k, i. e. nothing to update
	}

	//=====loop over all healthy residents=====//
	t_infec_k = househ[h].get_data("end_followup");
	for (std::size_t k : househ[h].susceptibles_list)
	{
		cumul_profile_lookuptab[from_to(i, k)] = compute_infec_profile_cumul(i, t_infec_k);
	}
}

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
double Model::get_infec_time(size_t i) const
{
	return indiv[i].get_augm("onset_time") - indiv[i].get_augm("incub_period");
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

    // compute infectivity profiles for each of the symptomatic and relative infectees
    for (size_t i = 0; i < num_indiv; i++)
    {
        store_augm_var_dependecies(i);
    }

    // find the first case for each household
    for (size_t h = 0; h < num_househ; h++)
    {
		// FIXME : for the moment there's a mess between set and assignemtn, it might not work because of this
        assign_end_followup(h);
        update_first_case(h);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////// Model : associated rates  ////////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the parameters controlling infectivity for each individual */
void Model::assign_rate_per_group_to_individual(size_t i)
{
	//=====rate 0 : Infectivity=====//
	if ((indiv[i].get_data("recov") == false) && (indiv[i].get_data("vacc_status") == 0))
	{ // unrecovered unvaccinated
		indiv[i].assign_rate_per_group("Infectivity", "i_naive");
	}
	else
	{ // vaccinated o recovered : Infectivity [REFERENCE]
		indiv[i].assign_rate_per_group("Infectivity", REFERENCE);
	}

	//=====rate 1 : Susceptiblity=====//
	if (indiv[i].get_data("vacc_status") == 0)
	{ // naive : Susceptiblity [REFERENCE]
		indiv[i].assign_rate_per_group("Susceptiblity", REFERENCE);
	}
	else
	{ // 1 dose
		indiv[i].assign_rate_per_group("Susceptiblity", "s_vacc");
	}
}

/* Assign the contact patterns between a couple of individuals. */
void Model::assign_rate_per_contact(std::size_t i, std::size_t k)
{
	bool i_is_child = (indiv[i].get_data("age") < get_value("child_age")); // children if younger than 12 y.o.
	bool i_is_adult = (i_is_child == false);							   // adults otherwise

	bool k_is_child = (indiv[k].get_data("age") < get_value("child_age")); // children if younger than 12 y.o.
	bool k_is_adult = (k_is_child == false);							   // adults otherwise

	if (i_is_child && k_is_child)
	{ // Child to Child
		contact_patterns[from_to(i, k)] = get_address("c_kid_to_kid");
	}
	else if (i_is_child && k_is_adult)
	{ // Child to Adult
		contact_patterns[from_to(i, k)] = get_address("c_kid_to_adu");
	}
	else if (i_is_adult && k_is_child)
	{ // Adult to Child
		contact_patterns[from_to(i, k)] = get_address("c_adu_to_kid");
	}
	else
	{ // Adult to Adult : Contact Pattern [REFERENCE]
		contact_patterns[from_to(i, k)] = get_reference();
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

/* Sum over different hazards. */
double Model::_sum_hazards_(double *profiles, std::size_t i)
{
	double output = 0.0;
	std::size_t h = indiv[i].home;

	// Risk of infection from other household members
	for (std::size_t k : househ[h].infectious_list)
	{
		if (k != i)
		{
			// 			contact pattern					infectivity			infectivity profile
			output += get_contact_from_to(k, i) * indiv[k].model_rate("Infectivity") * profiles[from_to(k, i)];
		}
	}

	return output * get_value("beta") * pow(get_value("median_size") / househ[h].house_size, get_value("delta"));
}

/* Compute the log-likelihood within household. */
double Model::compute_single_house_loglike(std::size_t h)
{
	double loglike_house = 0.0;
	double ti;
	double t0 = househ[h].get_data("start_followup");
	double tf = househ[h].get_data("end_followup");
	double comm_risk_inst, comm_risk_cumul;

	comm_risk_cumul = get_value("alpha") * (lin_idx_interp(tf, Cumul_EPI_Curve) - lin_idx_interp(t0, Cumul_EPI_Curve));

	// Contribution of healthy residents
	for (std::size_t i : househ[h].susceptibles_list)
	{
		loglike_house -= indiv[i].model_rate("Susceptibility") * (comm_risk_cumul + _sum_hazards_(cumul_profile_lookuptab, i));
	}

	// Contribution of symptomatic residents
	for (std::size_t i : househ[h].infectious_list)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : augm_with_model_set)
		{
			loglike_house += indiv[i].augm_model_log_pdf(augmID);
		}

		// Contribution of secondary cases
		if (i != househ[h].first_case)
		{
			ti = get_infec_time(i);
			comm_risk_inst = get_value("alpha") * lin_idx_interp(ti, Inst_EPI_Curve);
			loglike_house += log(indiv[i].model_rate("Susceptibility") * (comm_risk_inst + _sum_hazards_(inst_profile_lookuptab, i)));
			loglike_house -= indiv[i].model_rate("Susceptibility") * (comm_risk_cumul + _sum_hazards_(cumul_profile_lookuptab, i));
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

/* Poisson point process to simulate data. */
void Model::simulate_infection()
{
	double dt = 0.01;
	double proba, comm_risk_inst;

	// simulate the infection for each household
	for (std::size_t h = 0; h < num_househ; h++)
	{
		// fix all incub period at the beginning (all others should already be initialized)
		for (std::size_t i : househ[h].susceptibles_list)
		{
			indiv[i].init_augm("incub_period", rng);
		}

		// set all secondary cases to healthy
		std::vector<std::size_t> copy_list = househ[h].infectious_list;
		for (std::size_t i : copy_list)
		{
			if (i != househ[h].first_case)
			{
				indiv[i].set_status(HEALTHY);
				househ[h].change_status(i, HEALTHY);
			}
		}

		//=====simulate infection within household=====//
		for (double t = dt + househ[h].get_data("start_followup"); ((t < househ[h].get_data("end_followup")) && !househ[h].susceptibles_list.empty()); t += dt)
		{
			std::vector<std::size_t> copy_list = househ[h].susceptibles_list;
			for (std::size_t i : copy_list)
			{
				househ[h].change_status(i, SYMPTOMATIC);
				indiv[i].set_status(SYMPTOMATIC);

				indiv[i].set_augm("onset_time", t + uniform(rng) * dt + indiv[i].get_augm("incub_period"));
				store_augm_var_dependecies(i);
				comm_risk_inst = get_value("alpha") * lin_idx_interp(t, Inst_EPI_Curve);
				proba = 1 - exp(-indiv[i].model_rate("Susceptibility") * (comm_risk_inst + _sum_hazards_(inst_profile_lookuptab, i)) * dt);

				// refuse
				if (uniform(rng) >= proba)
				{
					househ[h].change_status(i, HEALTHY);
					indiv[i].set_status(HEALTHY);
				}
			}
		}
	}

	// Output formatter
	for (std::size_t i = 0; i < num_indiv; i++)
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