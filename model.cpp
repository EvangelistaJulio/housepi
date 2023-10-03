// model.cpp //

#include "model.h"

using namespace std;

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Comm_Infec;
extern std::vector<double> Inst_Comm_Infec;

//////////////////////////////////////////////////////////
// //////// Time Onset to Time of Transmission //////// //
//////////////////////////////////////////////////////////

/* The Time from Onset of Symptoms to Transmission (TOST) is modeled as a skew-logistic distribution (type I).
This distribution is also a function of the incubation time `t_incub`.
Note: symtoms can occur after infection (pre-symptomatic transmission), in other words `tost` can be negative.
----
Ref:
----
Ferretti, L. et al. The timing of COVID-19 transmission. 2020.09.04.20188516
Preprint at https://doi.org/10.1101/2020.09.04.20188516 (2020).*/

/* Probability distribution function */
double Onset_to_Transm::pdf(timetype tost, timetype t_incub)
{
	double out = 0.0;
	double a = get_arg("alpha");
	double m = get_arg("mu");
	double s = get_arg("sigma");
	double T = get_arg("tau");

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
	double a = get_arg("alpha");
	double m = get_arg("mu");
	double s = get_arg("sigma");
	double T = get_arg("tau");

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

/* Add individual */ // FIXME: base model
void Model::add_individual(Individual new_indiv)
{
	int status = new_indiv.get_status();

	if ((status != HEALTHY) && (status != SYMPTOMATIC))
	{
		std::cerr << "ERROR: status not in the model: " << status << endl;
		std::cerr << "Exit...\n"
				  << std::endl;
		std::exit(EXIT_FAILURE);
	}

	indiv.push_back(new_indiv);
	num_indiv++;
}

/* Assign community structure to individuals */
void Model::assign_community_structure()
{
	size_t househ_idx, size;
	std::unordered_map<int, Household> load_househ = {};

	//=====safe load households=====//
	househ.clear();
	for (size_t i = 0; i < num_indiv; i++)
	{
		househ_idx = indiv[i].get("HH_idx");
		// if househ_idx doesn't exist, append new household
		if (load_househ.find(househ_idx) == load_househ.end())
		{
			load_househ.emplace(househ_idx, Household());
			load_househ.at(househ_idx).set_index(househ_idx);
		}
		load_househ.at(househ_idx).add_individual(i, indiv[i].get_status());
	}
	for (const auto &pair : load_househ)
	{
		househ.push_back(pair.second);
	}
	num_househ = househ.size();
	load_househ.clear();

	//=====assign home to individuals=====//
	num_symptomatic = 0;
	for (size_t h = 0; h < num_househ; h++)
	{
		size = 0;
		for (size_t i : househ[h].symptomatic_list)
		{
			indiv[i].home = h;
			num_symptomatic += 1;
			size++;
		}
		for (size_t i : househ[h].healthy_list)
		{
			indiv[i].home = h;
			size++;
		}
		househ[h].size = size;
	}
}

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
		double init_value = 0.;

		std::string model_name = "gamma";
		kwargs_ext model_args = {
			{"shape", get_address("incub_shape")},
			{"scale", get_address("incub_scale")}};

		std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "normal";
		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", 0.5}};

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
		double init_value = (status == HEALTHY) ? indiv[i].get("study_period") : indiv[i].get("onset_time");

		std::string model_name = UNDEFINED;
		kwargs_ext model_args = {};

		std::string metropolis_strategy = (status == HEALTHY) ? UNDEFINED : "normal";

		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", get_value("onset_err")}};

		double metropolis_zero = (status == HEALTHY) ? INFINITY : init_value;
		kwargs range = {
			{"low", get_value("onset_low")},
			{"up", indiv[i].get("study_period") + 0.5}}; // FIXME onset_time == study_period!!!

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
	augm_to_be_inferred.push_back("onset_time");
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
		{"alpha", get_address("tost_alpha")},
		{"mu", get_address("tost_mu")},
		{"sigma", get_address("tost_sigma")},
		{"tau", get_address("tost_tau")}};

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
	size_t w;
	size_t h = indiv[i].home;
	timetype t_infec_i = get_infec_time(i);
	timetype t_infec_k;

	// Update cumulative and instantaneous person-to-person transmission rates
	// cumul_profile[idx + res * num_indiv] : from `idx` to `res`
	// cumul_profile[idx][res] : from `idx` to `res`

	//=====loop over all symptomatic residents=====//
	for (size_t k : househ[h].symptomatic_list)
	{
		t_infec_k = get_infec_time(k);

		if (t_infec_i < t_infec_k)
		{ // `i` : infector, `k` : infectee
			w = i + num_indiv * k;
			cumul_profile[w] = infectivity_cumul(i, t_infec_k);
			inst_profile[w] = infectivity_inst(i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			w = k + num_indiv * i;
			cumul_profile[w] = infectivity_cumul(k, t_infec_i);
			inst_profile[w] = infectivity_inst(k, t_infec_i);
		}
		// NOTE : else is the case idx == res
	}

	//=====loop over all healthy residents=====//
	// UPDATE : if study_period is equal for all residents, than this is constant for all
	for (size_t k : househ[h].healthy_list)
	{
		w = i + num_indiv * k;
		t_infec_k = get_infec_time(k); // == study_period
		cumul_profile[w] = infectivity_cumul(i, t_infec_k);
	}
}

/* Clear infectivity profiles */
void Model::clear_infectivity_profiles()
{
	cumul_profile.resize(pow(num_indiv, 2));
	inst_profile.resize(pow(num_indiv, 2));
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
	return get_value("alpha") * lin_idx_interp(t1, Inst_Comm_Infec);
}

/* Cumulative risk of infection from community */
double Model::community_risk_cumul(timetype t0, timetype t1)
{
	return get_value("alpha") * (lin_idx_interp(t1, Cumul_Comm_Infec) - lin_idx_interp(t0, Cumul_Comm_Infec));
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
}

/* Relative infectivity of symptomatic child cases versus symptomatic adult cases */
double Model::contact_modulation(size_t i, size_t k)
{
	double f = 1.0;

	// Child - Child
	if (indiv[i].is_child && indiv[k].is_child)
	{
		f *= get_value("rCC");
	} // Child - Adult
	else if ((indiv[i].is_child && indiv[k].is_adult) || (indiv[k].is_child && indiv[i].is_adult))
	{
		f *= get_value("rCA");
	}

	return f;
}

/* Model the susceptibility within household */
double Model::susceptibility(size_t i)
{
	double f = 1.0;

	if (indiv[i].is_vacc_0)
	{
		f *= get_value("rU"); // vaccination class 0
	}
	else if (indiv[i].is_vacc_1)
	{
		f *= get_value("r1V"); // vaccination class 1
	}
	else if (indiv[i].is_vacc_2)
	{
		f *= get_value("r2V"); // vaccination class 1
	}

	// Suceptibility of isolated people
	if (indiv[i].is_isol)
	{
		f *= get_value("rQ");
	}

	return f;
}

/* Model the transmissability within household */
double Model::transmissability(size_t i)
{
	size_t h = indiv[i].home;

	// Instantaneous per capita transmission rate
	double f = get_value("beta");

	// If transmission hazard depends on household size
	f *= pow(get_value("typical_size") / househ[h].size, get_value("delta"));

	return f;
}

/* beta risk */
double Model::risk_beta(std::vector<double> &rate_vector, size_t i)
{
	size_t w, h = indiv[i].home;
	double beta = 0.0;

	// Risk of infection from other household members
	for (size_t k : househ[h].symptomatic_list)
	{
		if (k != i)
		{
			w = k + i * num_indiv;
			beta += contact_modulation(k, i) * rate_vector[w];
		}
	}
	beta *= transmissability(i);
	return beta;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (likelihood)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Survival probability from time t0 of resident idx */
double Model::log_surv_proba(size_t i)
{
	size_t h = indiv[i].home;
	double beta = 0.0, output = 0.0;
	timetype t0 = househ[h].infec_time;
	timetype t_transm = get_infec_time(i);

	beta = risk_beta(cumul_profile, i);
	output = susceptibility(i) * (community_risk_cumul(t0, t_transm) + beta);
	output *= (-1.);

	return output;
}

/* logarithm of the infection probability */
double Model::log_infec_proba(size_t i)
{
	double beta = 0.0, output = 0.0;
	timetype t_transm = get_infec_time(i);

	beta = risk_beta(inst_profile, i);
	output = susceptibility(i) * (community_risk_inst(t_transm) + beta);
	output = std::log(1. - std::exp(-output));

	return output;
}

/* Compute the log-likelihood within household */
double Model::compute_loglike(size_t h)
{
	double loglike_house = 0.0;

	// Loop over all healthy residents
	for (size_t i : househ[h].healthy_list)
	{
		loglike_house += log_surv_proba(i);
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
			loglike_house += log_surv_proba(i); // FIXME : why not first_case?
			loglike_house += log_infec_proba(i);
		}
	}

	return loglike_house;
}