// model.cpp //

#include "model.h"

using namespace std;

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Comm_Infec;
extern std::vector<double> Inst_Comm_Infec;

//////////////////////////////////////////////////////////
// //////// Time Onset to Time of Transmission //////// //
//////////////////////////////////////////////////////////

/* Probability distribution function */
double Onset_to_Transm::pdf(timetype tost, timetype t_incub)
{
	double out = 0.0;
	double a = get_arg("tost_alpha"); // FIXME : gigantic bottleneck
	double m = get_arg("tost_mu");
	double s = get_arg("tost_sigma");
	double T = get_arg("tost_tau");

	double norm_const = (s / a) * (1. + (t_incub / T - 1.) * std::pow(1 + std::exp(m / s), -a) - t_incub / T * (std::pow(1. + std::exp((m + T) / s), -a)));

	if (tost < 0)
	{
		tost *= T / t_incub;
	}
	out = std::exp((m - tost) / s) * std::pow(1. + std::exp((m - tost) / s), -1. - a);

	return out / norm_const;
}

/* Cumulative distribution function */
double Onset_to_Transm::cdf(timetype tost, timetype t_incub)
{
	double out = 0.0;
	double a = get_arg("tost_alpha"); // FIXME : gigantic bottleneck
	double m = get_arg("tost_mu");
	double s = get_arg("tost_sigma");
	double T = get_arg("tost_tau");

	double norm_const = (s / a) * (1. + (t_incub / T - 1.) * std::pow(1 + std::exp(m / s), -a) - t_incub / T * (std::pow(1. + std::exp((m + T) / s), -a)));

	if (tost < 0)
	{
		// Integrate from -t_incub to tost
		out = t_incub / T * (std::pow(1. + std::exp((m - tost * T / t_incub) / s), -a) - std::pow(1. + std::exp((m + T) / s), -a));
	}
	else
	{
		// Integrate from -t_incub to 0
		out = t_incub / T * (std::pow(1. + std::exp(m / s), -a) - std::pow(1. + std::exp((m + T) / s), -a));
		// Integrate from 0 to tost
		out += std::pow(1. + std::exp((m - tost) / s), -a) - std::pow(1. + std::exp(m / s), -a);
	}

	return out / norm_const;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Model : base model /////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* First study date in the sample */
timetype Model::study_time_start()
{
	timetype output = 0;
	for (size_t i = 0; i < num_indiv; i++)
	{
		output = std::min(output, indiv[i].get("onset_time"));
	}

	return output;
}

/* Last study date in the sample */
timetype Model::study_time_end()
{
	timetype output = 0;
	for (size_t i = 0; i < num_indiv; i++)
	{
		output = std::max(output, indiv[i].get("study_period"));
	}

	return output;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// //////////////////////////////// Model : Augmented Variables //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Initialize augmented variables */
// FIXME: this needs to become more user friendly
void Model::_assign_augm_to_ind_(size_t i, std::mt19937_64 &rng)
{
	int status = get_status(i);

	//=====incub_period=====//
	{
		std::string augmID = "incub_period";
		// FLAG
		double init_value = (status == HEALTHY) ? 0.0 : indiv[i].get("incub_period");

		std::string model_name = "gamma";
		kwargs_ext model_args = {
			{"shape", get_address("incub_shape")},
			{"scale", get_address("incub_scale")}};

		// FLAG
		// std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "normal";
		std::string metropolis_strategy = UNDEFINED;

		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", 0.5}};

		// FLAG
		double metropolis_zero = INFINITY; // from old values

		kwargs range = {
			{"low", 0.0},
			{"up", get_value("max_incub_period")}};

		Augmented augm_var(
			init_value, model_name, model_args,
			metropolis_strategy, metropolis_args, metropolis_zero, range);

		// FLAG
		//(status == HEALTHY) ? indiv[i].add_augm(augmID, augm_var) : indiv[i].add_augm(augmID, augm_var, rng);
		indiv[i].add_augm(augmID, augm_var);
	}

	//=====onset_time=====//
	{
		std::string augmID = "onset_time";
		double init_value = (status == HEALTHY) ? indiv[i].get("study_period") : indiv[i].get("onset_time");

		std::string model_name = UNDEFINED;
		kwargs_ext model_args = {};

		std::string metropolis_strategy = UNDEFINED;
		kwargs metropolis_args = {};

		double metropolis_zero = INFINITY;
		kwargs range = {};

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
		_assign_augm_to_ind_(i, rng);
	}

	for (size_t i = 0; i < num_indiv; i++)
	{
		update_infectivity_profiles(i);
	}

	for (size_t h = 0; h < num_househ; h++)
	{
		update_first_case(h);
	}

	// FIXME (temporary) : this should be automatic
	augm_list.push_back("incub_period");
	augm_list.push_back("onset_time");
	augm_with_model.push_back("incub_period");
	// FLAG
	// augm_to_be_inferred.push_back("incub_period");
}

/* Time of infection */
timetype Model::get_infec_time(size_t i) const
{
	return indiv[i].get_augm("onset_time") - indiv[i].get_augm("incub_period");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////////// Model : Infectivity Profiles /////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

// FIXME : names overlapping...
void Model::initialize_infectivity_profile()
{
	kwargs_ext infec_prof_args = {
		{"tost_alpha", get_address("tost_alpha")},
		{"tost_mu", get_address("tost_mu")},
		{"tost_sigma", get_address("tost_sigma")},
		{"tost_tau", get_address("tost_tau")}};

	infectivity_profile = new Onset_to_Transm(infec_prof_args);
	clear_infectivity_profiles();
}

/* instanteneous infectivity profile */
double Model::infectivity_inst(size_t i, timetype t_transm)
{
	timetype tost = t_transm - indiv[i].get_augm("onset_time");
	timetype incub_period = indiv[i].get_augm("incub_period");
	return infectivity_profile->pdf(tost, incub_period);
}

/* cumulative infectivity profile */
double Model::infectivity_cumul(size_t i, timetype t_transm)
{
	timetype tost = t_transm - indiv[i].get_augm("onset_time");
	timetype incub_period = indiv[i].get_augm("incub_period");
	return infectivity_profile->cdf(tost, incub_period);
}

/* Update infectivity profiles given a resident idx */
void Model::update_infectivity_profiles(size_t i)
{
	size_t h = indiv[i].home;
	timetype t_infec_i = get_infec_time(i);
	timetype t_infec_k;

	// Update cumulative and instantaneous person-to-person transmission rates
	// NOTE : cumul_profile[idx][res] : from `idx` to `res`

	//=====loop over all symptomatic residents=====//
	for (size_t k : househ[h].symptomatic_list)
	{
		t_infec_k = get_infec_time(k);

		if (t_infec_i < t_infec_k)
		{ // `i` : infector, `k` : infectee
			inst_profile[i][k] = infectivity_inst(i, t_infec_k);
			cumul_profile[i][k] = infectivity_cumul(i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			inst_profile[k][i] = infectivity_inst(k, t_infec_i);
			cumul_profile[k][i] = infectivity_cumul(k, t_infec_i);
		}
		// NOTE : else is the case idx == res
	}

	//=====loop over all healthy residents=====//
	// UPDATE : if study_period is equal for all residents, than this is constant for all
	for (size_t k : househ[h].healthy_list)
	{
		t_infec_k = get_infec_time(k);						 // == study_period
		inst_profile[i][k] = infectivity_inst(i, t_infec_k); // FIXME : (temporary) not useful in
		cumul_profile[i][k] = infectivity_cumul(i, t_infec_k);
	}
}

/* Clear infectivity profiles */
void Model::clear_infectivity_profiles()
{
	cumul_profile.resize(num_indiv);
	inst_profile.resize(num_indiv);
	for (size_t j = 0; j < num_indiv; j++)
	{
		cumul_profile[j].resize(num_indiv);
		inst_profile[j].resize(num_indiv);
	}
}

/* Update household first case and first infection time */
void Model::update_first_case(size_t h)
{
	timetype t_infec;
	timetype new_t0 = INFINITY;
	size_t new_case0 = num_indiv; // NOTE: init to impossible values

	for (size_t i : househ[h].symptomatic_list)
	{
		t_infec = get_infec_time(i);
		new_t0 = std::min(t_infec, new_t0);
		if (t_infec == new_t0)
		{
			new_case0 = i;
		}
	}
	househ[h].infec_time = new_t0;
	househ[h].first_case = new_case0;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (modulators)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Instaneous risk of infection from community */
double Model::community_risk_inst(timetype t1)
{
	UNUSED(t1);
	return get_value("alpha");
}

/* Cumulative risk of infection from community */
double Model::community_risk_cumul(timetype t0, timetype t1)
{
	return get_value("alpha") * (t1 - t0);
}

/* Assign all the model classes starting from the data file */
void Model::pre_process_data()
{
	for (size_t i = 0; i < num_indiv; i++)
	{
		//=====Age class=====//
		indiv[i].is_child = (indiv[i].get("age") < 12);
		indiv[i].is_adult = (indiv[i].get("age") >= 12);

		//=====Vaccination status=====//
		indiv[i].is_isol = indiv[i].get("isol");
		indiv[i].is_vacc_0 = (indiv[i].get("recov")) || (indiv[i].get("vacc_status") == 1);
		indiv[i].is_vacc_1 = (indiv[i].get("vacc_status") == 2 && indiv[i].get("vacc_delay") < get_value("vacc_onset"));
		indiv[i].is_vacc_2 = (indiv[i].get("vacc_status") == 3 && indiv[i].get("vacc_delay") < get_value("vacc_onset"));
	}

	//=====clean up list of headers for output files=====//
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "index"));
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "status"));
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "onset_time"));
	// FLAG
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "incub_period"));
}

/* Relative infectivity of symptomatic child cases versus symptomatic adult cases */
double Model::contact_modulation(size_t i, size_t k)
{
	UNUSED(i);
	UNUSED(k);
	return 1.0;
}

/* Model the susceptibility within household */
double Model::susceptibility(size_t i)
{
	UNUSED(i);
	return 1.0;
}

/* Model the transmissability within household */
double Model::transmissability(size_t i)
{
	return get_value("beta") / househ[indiv[i].home].size;
}

/* beta risk */
double Model::risk_beta(std::vector<std::vector<double>> &rate_vector, size_t i)
{
	size_t h = indiv[i].home;
	double beta = 0.0;

	// Risk of infection from other household members
	for (size_t k : househ[h].symptomatic_list)
	{
		if (k != i)
		{
			beta += contact_modulation(k, i) * rate_vector[k][i];
		}
	}
	beta *= transmissability(i);
	return beta;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Survival probability from time t0 of resident idx */
double Model::risk_surv(size_t i)
{
	size_t h = indiv[i].home;
	timetype t0 = househ[h].infec_time;
	timetype t_transm = get_infec_time(i);
	return susceptibility(i) * (community_risk_cumul(t0, t_transm) + risk_beta(cumul_profile, i));
}

/* logarithm of the infection probability */
double Model::risk_infec(size_t i)
{
	timetype t_transm = get_infec_time(i);
	return susceptibility(i) * (community_risk_inst(t_transm) + risk_beta(inst_profile, i));
}

/* Compute the log-likelihood within household */
double Model::compute_loglike(size_t h)
{
	double loglike_house = 0.0;

	// Loop over all healthy residents
	for (size_t i : househ[h].healthy_list)
	{
		loglike_house -= risk_surv(i);
	}

	// Loop over all symptomatic residents
	for (size_t i : househ[h].symptomatic_list)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : augm_with_model)
		{
			loglike_house += indiv[i].augm_model_log_pdf(augmID);
		}

		// Contribution of secondary cases
		if (i != househ[h].first_case)
		{
			loglike_house -= risk_surv(i); // NOTE : risk_surv(first_case) == 0
			loglike_house += std::log(1. - std::exp(-risk_infec(i)));
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
	timetype end_followup;
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	double proba;

	// simulate the infection for each household
	for (size_t h = 0; h < num_househ; h++)
	{
		// fix all incub period at the beginning
		// all others should already be initialized
		for (size_t i : househ[h].healthy_list)
		{
			indiv[i].init_augm("incub_period", rng);
		}

		// find maximum study periodo within the household
		end_followup = 0;
		std::vector<size_t> copy_list = househ[h].symptomatic_list;
		for (size_t i : copy_list)
		{
			// set all secondary cases to healthy
			if (i != househ[h].first_case)
			{
				indiv[i].set_status(HEALTHY);
				househ[h].change_status(i, HEALTHY);
			}
			end_followup = std::max(end_followup, indiv[i].get("study_period"));
		}

		//=====simulate infection within household=====//
		for (timetype t = househ[h].infec_time; ((t < end_followup) && !househ[h].healthy_list.empty()); t += dt)
		{
			std::vector<size_t> copy_list = househ[h].healthy_list;
			for (size_t i : copy_list)
			{
				// CHECK : it should just be at time t to move to t + dt (Euler's scheme)
				indiv[i].set_augm("onset_time", t + dt + indiv[i].get_augm("incub_period")); // unif(rng) *
				update_infectivity_profiles(i);

				// proba = (1. - std::exp(-risk_infec(i))) * dt; // approximate with mean value theorem
				proba = 1. - std::exp(-risk_infec(i) * dt); // approximate with Paolo's expertise

				// FLAG
				// double u = unif(rng);
				// cerr << indiv[i].get_index() << " " << t << " " << u << " " << proba << endl;
				// if (u < proba)
				// FLAG

				if (unif(rng) < proba)
				{
					househ[h].change_status(i, SYMPTOMATIC);
					indiv[i].set_status(SYMPTOMATIC);
				}
			}

			for (size_t i : househ[h].symptomatic_list)
			{
				update_infectivity_profiles(i);
			}
		}
	}

	// formatter
	for (size_t i = 0; i < num_indiv; i++)
	{
		if (indiv[i].get_status() == SYMPTOMATIC)
		{
			if (indiv[i].get_augm("onset_time") > indiv[i].get("study_period"))
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
			// FLAG
			// indiv[i].set_augm("incub_period", 0);
		}
	}
}