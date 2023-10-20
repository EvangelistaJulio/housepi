// model.cpp //

#include "model.h"

using namespace std;

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Comm_Infec;
extern std::vector<double> Inst_Comm_Infec;


///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// Individual : Methods for Model ////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Time of infection */
timetype Individual::get_infec_time() const
{
	return get_augm("onset_time") - get_augm("incub_period");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Time Onset to Time of Transmission //////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* The Time from Onset of Symptoms to Transmission (TOST) is modeled as a skew-logistic distribution (type I).
This distribution is also a function of the incubation time `t_incub`.
Note: symtoms can occur after infection (pre-symptomatic transmission), in other words `tost` can be negative.
----
Ref:
----
Ferretti, L. et al. The timing of COVID-19 transmission. 2020.09.04.20188516
Preprint at https://doi.org/10.1101/2020.09.04.20188516 (2020).*/

/* Probability distribution function for the Infectivity Profile. */
double Onset_to_Transm::pdf(timetype tost, timetype t_incub)
{
	double a = get_arg(0);
	double m = get_arg(1);
	double s = get_arg(2);
	double T = get_arg(3);

	double norm_pdf = (s / a) * (1. + (t_incub / T - 1.) * std::pow(1 + std::exp(m / s), -a) - (t_incub / T) * (std::pow(1. + std::exp((m + T) / s), -a)));

	if (tost < 0)
	{
		tost *= (T / t_incub);
	}
	return std::exp((m - tost) / s) * std::pow(1. + std::exp((m - tost) / s), -1. - a) / norm_pdf;
}

/* Cumulative distribution function forthe Infectivity Profile. */
double Onset_to_Transm::cdf(timetype tost, timetype t_incub)
{
	double a = get_arg(0);
	double m = get_arg(1);
	double s = get_arg(2);
	double T = get_arg(3);

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

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////// Model : Augmented Variables //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize augmented variables */
// FIXME: this needs to become more user friendly
void Model::_assign_augm_to_individual_(size_t i, std::mt19937_64 &rng)
{
	int status = get_status(i);

	//=====incub_period=====//
	{
		std::string augmID = "incub_period";
		double init_value = 0.;

		std::string model_name = "gamma";
		kwargs_ext model_args = {
			{"shape", get_address("incub_shape")},
			{"scale", get_address("incub_scale")}};

		// FIXME : verify that metropolis_strategy = prior
		/*
		std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "normal";
		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", 0.5}};
		*/
		std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "prior";
		kwargs metropolis_args = {};

		double metropolis_zero = INFINITY; // from old value
		kwargs range = {
			{"low", 0.0},
			{"up", get_value("max_incub_period")}};

		Augmented augm_var(
			init_value, model_name, model_args,
			metropolis_strategy, metropolis_args, metropolis_zero, range);
		(status == HEALTHY) ? indiv[i].add_augm(augmID, augm_var) : indiv[i].add_augm(augmID, augm_var, rng);
	}

	//=====onset_time=====//
	{
		std::string augmID = "onset_time";
		double init_value = (status == HEALTHY) ? indiv[i].get_data("study_period") : indiv[i].get_data("onset_time");

		std::string model_name = UNDEFINED;
		kwargs_ext model_args = {};

		std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "normal";

		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", get_value("onset_err")}};

		double metropolis_zero = (status == HEALTHY) ? INFINITY : init_value;
		kwargs range = {
			{"low", get_value("onset_low")},
			{"up", indiv[i].get_data("study_period") + 0.1}}; // FIXME onset_time == study_period!!!

		Augmented augm_var(
			init_value, model_name, model_args,
			metropolis_strategy, metropolis_args, metropolis_zero, range);
		indiv[i].add_augm(augmID, augm_var);
	}
}

// Initialize augmented variables and infectivity profiles
void Model::assign_augm_data(std::mt19937_64 &rng)
{
	// init augmented data for all residents (user defined)
	for (size_t i = 0; i < num_indiv; i++)
	{
		_assign_augm_to_individual_(i, rng);
	}

	// FIXME (temporary) : this should be automatic
	augm_list.push_back("incub_period");
	augm_list.push_back("onset_time");
	augm_to_be_inferred.push_back("incub_period");
	augm_to_be_inferred.push_back("onset_time");
	augm_with_model.push_back("incub_period");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Model : Infectivity Profiles /////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize the infectivity profiles. */
void Model::initialize_infec_profiles()
{
	kwargs_ext infec_prof_args = {
		{"tost_alpha", get_address("tost_alpha")},
		{"tost_mu", get_address("tost_mu")},
		{"tost_sigma", get_address("tost_sigma")},
		{"tost_tau", get_address("tost_tau")}};

	infec_profile_func = new Onset_to_Transm(infec_prof_args);

	// produce the empty matrix for the infectivity profiles
	inst_profile = new double[num_indiv * num_indiv];
	cumul_profile = new double[num_indiv * num_indiv];

	// compute infectivity profiles for each of the symptomatic and relative infectees
	for (size_t i = 0; i < num_indiv; i++)
	{
		update_infec_profile(i);
	}
}

/* Instanteneous infectivity profile. */
double Model::infec_profile_inst(size_t i, timetype t_transm)
{
	timetype tost = t_transm - indiv[i].get_augm("onset_time");
	timetype incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func->pdf(tost, incub_period);
}

/* Cumulative infectivity profile. */
double Model::infec_profile_cumul(size_t i, timetype t_transm)
{
	timetype tost = t_transm - indiv[i].get_augm("onset_time");
	timetype incub_period = indiv[i].get_augm("incub_period");
	return infec_profile_func->cdf(tost, incub_period);
}

/* Update cumulative and instantaneous person-to-person infectivity profiles according to an infector. */
void Model::update_infec_profile(size_t i)
{
	size_t h = indiv[i].home;
	timetype t_infec_i = indiv[i].get_infec_time();
	timetype t_infec_k;

	//=====loop over all symptomatic residents=====//
	for (size_t k : househ[h].infectors_list)
	{
		t_infec_k = indiv[k].get_infec_time();

		if (t_infec_i < t_infec_k)
		{ // `i` : infector, `k` : infectee
			inst_profile[from_to(i, k)] = infec_profile_inst(i, t_infec_k);
			cumul_profile[from_to(i, k)] = infec_profile_cumul(i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			inst_profile[from_to(k, i)] = infec_profile_inst(k, t_infec_i);
			cumul_profile[from_to(k, i)] = infec_profile_cumul(k, t_infec_i);
		}
		// NOTE : else is the case idx == res
	}

	//=====loop over all healthy residents=====//
	for (size_t k : househ[h].susceptibles_list)
	{
		t_infec_k = indiv[k].get_infec_time(); // == study_period
		cumul_profile[from_to(i, k)] = infec_profile_cumul(i, t_infec_k);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (modulators)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Instaneous risk of infection from community */
double Model::community_risk_inst(timetype t1)
{
	return get_value("alpha") * lin_idx_interp(t1, Inst_Comm_Infec);
}

/* Cumulative risk of infection from community */
double Model::community_risk_cumul(timetype t0, timetype t1)
{
	return get_value("alpha") * (lin_idx_interp(t1, Cumul_Comm_Infec) - lin_idx_interp(t0, Cumul_Comm_Infec));
}

/* Assign the parameters controlling infectivity and susceptibility for each individual */
// UPDATE : how to make this user-friendly
void Model::pre_process_data()
{
	const double *address;
	bool i_is_child, k_is_child;
	bool i_is_adult, k_is_adult;

	for (size_t i = 0; i < num_indiv; i++)
	{
		i_is_child = (indiv[i].get_data("age") < 12);
		i_is_adult = (i_is_child == false) && (indiv[i].get_data("sex") == 1);

		// ----- Susceptiblity ----- //
		if (indiv[i].get_data("isol") == true)
		{ // isolated
			indiv[i].set_suscept(get_address("s_isol"));
		}
		else
		{ // unisolated
			if (indiv[i].get_data("recov") == true)
			{ // recovered
				indiv[i].set_suscept(get_address("s_recov"));
			}
			else
			{ // unrecovered
				if (indiv[i].get_data("vacc_status") == 0)
				{ // unvaccinated : Susceptibility [REFERENCE]
					indiv[i].set_suscept(get_reference());
				}
				else if (indiv[i].get_data("vacc_status") == 1)
				{ // unisolated 1 dose
					indiv[i].set_suscept(get_address("s_vacc_1"));
				}
				else if ((indiv[i].get_data("vacc_status") == 2))
				{ // unisolated 2 doses - less vaccination onset
					indiv[i].set_suscept(get_address(("s_vacc_2")));
				}
				else
				{ // unisolated 3 doses
					indiv[i].set_suscept(get_address(("s_vacc_3")));
				}
			}
		}

		// ----- Infectivity ----- //
		if ((indiv[i].get_data("recov") == false) && (indiv[i].get_data("vacc_status") == 0))
		{ // unrecovered unvaccinated
			indiv[i].set_infecti(get_address("i_naive"));
		}
		else
		{ // vaccinated o recovered : Infectivity (REFERENCE)
			indiv[i].set_infecti(get_reference());
		}

		// ===== Contact Patterns ==== //

		for (size_t k = 0; k < i; k++)
		{
			k_is_child = (indiv[i].get_data("age") < 12);
			k_is_adult = (k_is_child == false);

			if (i_is_child && k_is_child)
			{ // Child - Child
				address = get_address("c_child_child");
			}
			else if ((i_is_child && k_is_adult) || (k_is_child && i_is_adult))
			{ // Child - Adult
				address = get_address("c_child_adult");
			}
			else
			{ // Adult - Adult : Contact Pattern(REFERENCE)
				address = get_reference();
			}
			contact_patterns[from_to(i, k)] = address;
			contact_patterns[from_to(k, i)] = address;
		}
	}

	//=====clean up list of headers for output files=====//
	// UPDATE : this shall be automatic
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "index"));
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "status"));
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "onset_time"));
}

/* Model the transmissability. */
double Model::transmissability(size_t h)
{
	return get_value("beta") * pow(get_value("typical_size") / househ[h].size, get_value("delta"));
}

/* Sum over different hazards. */
double Model::_sum_hazards_(double *profiles, size_t i)
{
	double output = 0.0;
	size_t h = indiv[i].home;

	// Risk of infection from other household members
	for (size_t k : househ[h].infectors_list)
	{
		if (k != i)
		{
			output += (*contact_patterns[from_to(k, i)]) * indiv[k].get_infecti() * profiles[from_to(k, i)];
		}
	}

	return output * transmissability(h) ;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Risk of infection. */
double Model::risk_infec(size_t i)
{
	timetype t_transm = indiv[i].get_infec_time();
	return indiv[i].get_suscept() * (community_risk_inst(t_transm) + _sum_hazards_(inst_profile, i));
}

/* Risk of surviving. */
double Model::risk_surv(size_t i)
{
	timetype t0 = househ[indiv[i].home].infec_time;
	timetype t_transm = indiv[i].get_infec_time();
	return indiv[i].get_suscept() * (community_risk_cumul(t0, t_transm) + _sum_hazards_(cumul_profile, i));
}

/* Compute the log-likelihood within household */
double Model::compute_loglike(size_t h)
{
	double loglike_house = 0.0;

	// Loop over all healthy residents
	for (size_t i : househ[h].susceptibles_list)
	{
		loglike_house -= risk_surv(i);
	}

	// Loop over all symptomatic residents
	for (size_t i : househ[h].infectors_list)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : augm_with_model)
		{
			loglike_house += indiv[i].augm_model_log_pdf(augmID);
		}

		// Contribution of secondary cases
		if (i != househ[h].first_case)
		{
			loglike_house += std::log(risk_infec(i)) - risk_surv(i);
		}
	}

	return loglike_house;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (simulation)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate data. */
void Model::simulate_infection(timetype dt, std::mt19937_64 &rng)
{
	double proba;
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

		// find maximum study periodo within the household
		end_followup = 0;
		std::vector<size_t> copy_list = househ[h].infectors_list;
		for (size_t i : copy_list)
		{
			// set all secondary cases to healthy
			if (i != househ[h].first_case)
			{
				indiv[i].set_status(HEALTHY);
				househ[h].change_status(i, HEALTHY);
			}
			end_followup = std::max(end_followup, indiv[i].get_data("study_period"));
		}

		//=====simulate infection within household=====//
		for (timetype t = dt + househ[h].infec_time; ((t < end_followup) && !househ[h].susceptibles_list.empty()); t += dt)
		{
			std::vector<size_t> copy_list = househ[h].susceptibles_list;
			for (size_t i : copy_list)
			{
				househ[h].change_status(i, SYMPTOMATIC);
				indiv[i].set_status(SYMPTOMATIC);

				indiv[i].set_augm("onset_time", t + unif(rng) * dt + indiv[i].get_augm("incub_period"));
				update_infec_profile(i);
				proba = 1 - exp(-risk_infec(i) * dt);

				// cerr << indiv[i].get_index() << " " << t << " " << proba << endl;

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