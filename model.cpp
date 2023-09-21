// model.cpp //

#include "model.hpp"

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

//////////////////////////////////////////////////////////
// /////////////////// User Model ///////////////////// //
//////////////////////////////////////////////////////////

// Add individual
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

	// safe load households
	for (const auto &pair : load_househ)
	{
		househ.push_back(pair.second);
	}
	num_househ = househ.size();

	// assign house to individuals
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
	load_househ.clear();
}

/* First study date in the sample */
timetype Model::first_study_date()
{
	timetype output = 0;
	for (size_t i = 0; i < num_indiv; i++)
	{
		output = std::min(output, indiv[i].get("onset_time"));
	}

	return output;
}

/* Last study date in the sample */
timetype Model::last_study_date()
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
// UPDATE: this needs to become more user friendly
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

		double value_zero = INFINITY; // from old value
		kwargs range = {{"low", 0.0}, {"up", get_value("max_incub_period")}};

		Augmented augm_var(init_value, model_name, model_args, metropolis_strategy, metropolis_args, value_zero, range);
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

		double value_zero = init_value; // from init value
		kwargs range = {{"low", get_value("onset_low")}, {"up", indiv[i].get("study_period")}};

		Augmented augm_var(init_value, model_name, model_args, metropolis_strategy, metropolis_args, value_zero, range);
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

void Model::assign_infec_profile()
{
	kwargs_ext infec_prof_args = {
		{"alpha", get_address("tost_alpha")},
		{"mu", get_address("tost_mu")},
		{"sigma", get_address("tost_sigma")},
		{"tau", get_address("tost_tau")}};

	infectivity_profile = new Onset_to_Transm(infec_prof_args);
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
			cumul_profile[i][k] = infectivity_cumul(i, t_infec_k);
			inst_profile[i][k] = infectivity_inst(i, t_infec_k);
		}
		else if (t_infec_k < t_infec_i)
		{ // `k` : infector, `i` : infectee
			cumul_profile[k][i] = infectivity_cumul(k, t_infec_i);
			inst_profile[k][i] = infectivity_inst(k, t_infec_i);
		}
		// NOTE : else is the case idx == res
	}

	//=====loop over all healthy residents=====//
	// UPDATE : if study_period is equal for all residents, than this is constant for all
	for (size_t k : househ[h].healthy_list)
	{
		t_infec_k = get_infec_time(k); // == study_period //FIXME : check
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

/* Update household first case */
void Model::update_first_case(size_t h)
{
	timetype t_infec, new_t0 = INFINITY;
	size_t new_case0 = househ[h].size; // NOTE: init to impossible values

	for (size_t i : househ[h].symptomatic_list)
	{
		// Find smaller infection time
		t_infec = get_infec_time(i);
		new_t0 = std::min(t_infec, new_t0);

		// Identify first case
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
		indiv[i].is_child = (indiv[i].get("age") == 0);
		indiv[i].is_mother = (indiv[i].get("age") == 1 && indiv[i].get("sex") == 1);
		indiv[i].is_father = (indiv[i].get("age") == 1 && indiv[i].get("sex") == 2);
		indiv[i].is_recov = (indiv[i].get("recov") && (indiv[i].get("vacc_status") == 1 || indiv[i].get("vacc_status") == 2));
		indiv[i].is_isol = indiv[i].get("isol");
		indiv[i].is_vacc_0 = (indiv[i].get("vacc_status") == 1);
		indiv[i].is_vacc_1 = (indiv[i].get("vacc_status") == 2 && indiv[i].get("vacc_delay") < get_value("vacc_onset"));
		indiv[i].is_vacc_1_V_onset = (indiv[i].get("vacc_status") == 2 && indiv[i].get("vacc_delay") >= get_value("vacc_onset"));
		indiv[i].is_vacc_2 = (indiv[i].get("vacc_status") == 3 && indiv[i].get("vacc_delay") < get_value("vacc_onset"));
	}
	// clean up the list of headers for output files
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "index"));
	data_headers.erase(std::find(data_headers.begin(), data_headers.end(), "status"));
}

// Relative infectivity of symptomatic child cases versus symptomatic adult cases
double Model::contact_modulation(size_t i, size_t k)
{
	double f = 1.0;

	// Child (2) - Child (2)
	if (indiv[i].is_child && indiv[k].is_child)
	{
		f *= get_value("rCC");

	} // Child (2) - Mother (0)
	else if ((indiv[i].is_child && indiv[k].is_mother) || (indiv[k].is_child && indiv[i].is_mother))
	{
		f *= get_value("rMC");

	} // Child (2) - Father (1)
	else if ((indiv[i].is_child && indiv[k].is_father) || (indiv[k].is_child && indiv[i].is_father))
	{
		f *= get_value("rFC");
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
	else if (indiv[i].is_vacc_1_V_onset)
	{
		f *= get_value("r1V90"); // vaccination class 1V90
	}
	else if (indiv[i].is_vacc_2)
	{
		f *= get_value("r2V"); // vaccination class 1
	}

	// Suceptibility of children
	if (indiv[i].is_child)
	{
		f *= get_value("rC");
	}

	// Suceptibility of isolated people
	if (indiv[i].is_isol)
	{
		f *= get_value("rQ");
	}

	// Suceptibility of recovered people
	if (indiv[i].is_recov)
	{
		f *= get_value("rR");
	}

	return f;
}

/* Model the transmissability within household */
double Model::transmissability(size_t i)
{
	double f = 1.0;
	size_t size = househ[indiv[i].home].size;

	// Instantaneous per capita transmission rate
	f *= get_value("beta");

	// If transmission hazard depends on household size
	if (get_infer("delta"))
	{
		f /= pow(size / get_value("typical_size"), get_value("delta"));
	}

	return f;
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
	for (size_t i = 0; i < househ[h].size; i++)
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
			loglike_house += log_infec_proba(i);
		}
	}

	return loglike_house;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////// Model : Methods (simulation)  /////////////////////////////// //
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Poisson point process to simulate data */
void Model::simulate_data(timetype delta_t, std::mt19937_64 &rng)
{
	size_t h;
	std::uniform_real_distribution<double> unif(0.0, 1.0);

	// initialize all individuals to healthy
	for (size_t i = 0; i < num_indiv; i++)
	{
		if (indiv[i].get_status() != HEALTHY)
		{
			househ[indiv[i].home].change_status(i, HEALTHY);
			indiv[i].set_status(HEALTHY);
		}

		indiv[i].init_augm("incub_period", rng); // UPDATE : merge with rand init
		indiv[i].set_augm("onset_time", INFINITY);
	}

	// initialize households
	for (size_t h = 0; h < num_househ; h++)
	{
		househ[h].infec_time = INFINITY;
		househ[h].first_case = househ[h].size;
	}
	clear_infectivity_profiles();

	timetype first_infec_day = INFINITY;
	for (timetype infec_time = first_study_date(); infec_time <= last_study_date(); infec_time += delta_t)
	{
		for (size_t i = 0; i < num_indiv; i++)
		{
			// UPDATE : state check improvable
			if (indiv[i].get_status() == HEALTHY)
			{
				// UPDATE : buffer-like construction could be helpful
				// NOTE : assign infec_time from onset_time
				indiv[i].set_augm("onset_time", infec_time + indiv[i].get_augm("incub_period")); 

				if (unif(rng) < std::exp(log_infec_proba(i)))
				{
					indiv[i].set_status(SYMPTOMATIC);
					h = indiv[i].home;
					househ[h].change_status(i, SYMPTOMATIC);
					update_first_case(h);
					update_infectivity_profiles(i);
					first_infec_day = std::min(infec_time, first_infec_day); // FIXME : sub-optimal
				}
				else
				{
					indiv[i].set_augm("onset_time", INFINITY);
				}
			}
		}
	}
}