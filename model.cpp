// model.cpp //

#include "model.hpp"

using namespace std;

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Community_Infec;
extern std::vector<double> Inst_Community_Infec;

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
		tost = tost * T / t_incub;
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

/* Initialize augmented variables */
// UPDATE: this needs to become more user friendly
void Model::_assign_augm_to_ind_(Household *H, size_t idx, std::mt19937_64 &rng)
{
	int status = H->get_status(idx);

	// --- incub_period --- //
	{
		std::string augmID = "incub_period";
		double init_value = 0.;

		std::string model_name = (status == HEALTHY) ? UNDEF : "gamma";
		kwargs_ext model_args = {
			{"shape", get_address("incub_shape")},
			{"scale", get_address("incub_scale")}};

		std::string metropolis_strategy = (status == HEALTHY) ? UNDEF : "normal";
		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", 0.5}};

		double value_zero = INFINITY; // from old value
		kwargs range = {{"low", 0.0}, {"up", get_value("max_incub_period")}};

		H->h_residents[idx].add_augm(augmID, init_value, model_name, model_args, metropolis_strategy, metropolis_args, value_zero, range);

		// FIXME : this should be automatic
		if (status != HEALTHY)
		{
			H->h_residents[idx].init_augm(augmID, rng);
		}
	}

	// --- onset_time --- //
	{
		std::string augmID = "onset_time";
		double init_value = (status == HEALTHY) ? H->h_residents[idx].i_study_period : H->h_residents[idx].i_onset_time;

		std::string model_name = UNDEF;
		kwargs_ext model_args = {};

		std::string metropolis_strategy = (status == HEALTHY) ? UNDEF : "normal";

		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", get_value("onset_err")}};

		double value_zero = init_value; // from init value
		kwargs range = {{"low", get_value("onset_low")}, {"up", H->h_residents[idx].i_study_period}};

		H->h_residents[idx].add_augm(augmID, init_value, model_name, model_args, metropolis_strategy, metropolis_args, value_zero, range);
	}
}

// Initialize augmented variables and infectivity profiles
void Model::assign_augm_data(Household *H, std::mt19937_64 &rng)
{
	// init augmented data for all residents (user defined)
	for (size_t res = 0; res < H->h_size; res++)
	{
		_assign_augm_to_ind_(H, res, rng);
	}
	update_first_case(H);

	for (size_t res = 0; res < H->h_size; res++)
	{
		update_infectivity_profiles(H, res);
	}
}

/* Time of infection */
timetype Model::get_infec_time(Household *H, size_t idx) const
{
	return H->get_augm("onset_time", idx) - H->get_augm("incub_period", idx);
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
double Model::infectivity_inst(Household *H, size_t idx, timetype t_transm)
{
	timetype tost = t_transm - H->get_augm("onset_time", idx);
	timetype incub_period = H->get_augm("incub_period", idx);
	return infectivity_profile->pdf(tost, incub_period);
}

/* cumulative infectivity profile */
double Model::infectivity_cumul(Household *H, size_t idx, timetype t_transm)
{
	timetype tost = t_transm - H->get_augm("onset_time", idx);
	timetype incub_period = H->get_augm("incub_period", idx);
	return infectivity_profile->cdf(tost, incub_period);
}

/* Update infectivity profiles given a resident idx */
void Model::update_infectivity_profiles(Household *H, size_t idx)
{
	timetype t_infec_idx = get_infec_time(H, idx);
	timetype t_infec_res;

	// Update cumulative and instantaneous person-to-person transmission rates
	// NOTE : h_cumul_profile[idx][res] : from `idx` to `res`

	// loop over all non healthy residents
	for (size_t res : H->h_symptomatic_list)
	{
		t_infec_res = get_infec_time(H, res);
		if (t_infec_idx < t_infec_res)
		{ // `idx` : infector, `res` : infectee
			H->h_cumul_profile[idx][res] = infectivity_cumul(H, idx, t_infec_res);
			H->h_inst_profile[idx][res] = infectivity_inst(H, idx, t_infec_res);
		}
		else if (t_infec_res < t_infec_idx)
		{ // `res` : infector, `idx` : infectee
			H->h_cumul_profile[res][idx] = infectivity_cumul(H, res, t_infec_idx);
			H->h_inst_profile[res][idx] = infectivity_inst(H, res, t_infec_idx);
		}
		// NOTE : else is the case idx == res
	}

	// loop over all healthy residents
	// UPDATE : if study_period is equal for all residents, than this is constant for all
	for (size_t res : H->h_healthy_list)
	{
		t_infec_res = get_infec_time(H, res); // == study_period
		H->h_cumul_profile[idx][res] = infectivity_cumul(H, idx, t_infec_res);
	}
}

/* Initialize infectivity profiles */
void Model::init_infectivity_profiles(Household *H)
{
	// initialize 2-dim vectors h_cumul_profile and h_inst_profile
	H->h_cumul_profile.resize(H->h_size);
	H->h_inst_profile.resize(H->h_size);

	// loop over all residents
	for (size_t i = 0; i < H->h_size; i++)
	{
		H->h_cumul_profile[i].resize(H->h_size);
		H->h_inst_profile[i].resize(H->h_size);
	}
}

/* Update household first case */
void Model::update_first_case(Household *H)
{
	timetype t_infec, new_t0 = INFINITY;
	size_t new_case0 = H->h_size; // NOTE: init to impossible values

	for (size_t i : H->h_symptomatic_list)
	{
		// Find smaller infection time
		t_infec = get_infec_time(H, i);
		new_t0 = std::min(t_infec, new_t0);

		// Identify first case
		if (t_infec == new_t0)
		{
			new_case0 = i;
		}
	}

	H->h_infec_time = new_t0;
	H->h_first_case = new_case0;
}

/////////////////////////////////////////////////////
//  Household : Methods (likelihood)
/////////////////////////////////////////////////////

double Model::risk_beta(Household *H, const std::vector<std::vector<double>> &rate_vector, size_t idx)
{
	double beta_i = 0.0, beta = 0.0;

	// Risk of infection from other household members
	for (size_t res : H->h_symptomatic_list)
	{
		if (idx != res)
		{
			beta_i = rate_vector[res][idx];
			beta_i *= contact_modulation(H, res, idx);
			beta += beta_i;
		}
	}
	beta *= transmissability(H);
	return beta;
}

/* Survival probability from time t0 of resident idx */
double Model::log_surv_proba(Household *H, size_t idx)
{
	double beta = 0.0, output = 0.0;
	timetype t0 = H->h_infec_time;
	timetype t_transm = get_infec_time(H, idx);

	beta = risk_beta(H, H->h_cumul_profile, idx);
	output = susceptibility(H, idx) * (community_risk_cumul(t0, t_transm) + beta);
	output *= (-1.);

	return output;
}

/* logarithm of the infection probability */
double Model::log_infec_proba(Household *H, size_t idx)
{
	double beta = 0.0, output = 0.0;
	timetype t_transm = get_infec_time(H, idx);

	beta = risk_beta(H, H->h_inst_profile, idx);
	output = susceptibility(H, idx) * (community_risk_inst(t_transm) + beta);
	output = std::log(1. - std::exp(-output));

	return output;
}

/* Compute the log-likelihood within household */
double Model::compute_loglike(Household *H)
{
	double loglike_house = 0.0;

	// Loop over all healthy residents
	// NOTE: no augmented variables
	for (size_t res = 0; res < H->h_size; res++)
	{
		loglike_house += log_surv_proba(H, res);

		// FLAG
		//std::cerr << "surv (" << res << ") = " << log_surv_proba(H, res) << std::endl;
		// FLAG
	}

	// Loop over all symptomatic residents
	for (size_t res : H->h_symptomatic_list)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : H->get_augm_with_model(res))
		{
			loglike_house += H->augm_model_log_pdf(augmID, res);

			// FLAG
			//std::cerr << "augm (" << res << ") = " << H->augm_model_log_pdf(augmID, res) << std::endl;
			// FLAG
		}

		// Contribution of secondary cases
		if (res != H->h_first_case)
		{
			loglike_house += log_infec_proba(H, res);

			// FLAG
			//std::cerr << "infec (" << res << ") = " << log_infec_proba(H, res) << std::endl;
			// FLAG
		}
	}

	// FLAG
	//std::exit(EXIT_FAILURE);
	// FLAG

	return loglike_house;
}

/////////////////////////////////////////////////////
///////////  Household : Methods (modulators) ///////
/////////////////////////////////////////////////////

/* Instaneous risk of infection from community */
double Model::community_risk_inst(timetype t1)
{
	double output = 0.0;
	output = lin_idx_interp(t1, Inst_Community_Infec);
	output *= get_value("alpha");

	return output;
}

/* Cumulative risk of infection from community */
double Model::community_risk_cumul(timetype t0, timetype t1)
{
	double output = 0.0;
	output = lin_idx_interp(t1, Cumul_Community_Infec) - lin_idx_interp(t0, Cumul_Community_Infec);
	output *= get_value("alpha");

	return output;
}

// Relative infectivity of symptomatic child cases versus symptomatic adult cases
double Model::contact_modulation(Household *H, size_t infectee, size_t idx)
{
	double factor = 1.0;
	size_t i = infectee;

	// Child (2) - Child (2)
	if (get_infer("rCC") && H->h_residents[idx].i_age == 0 && H->h_residents[i].i_age == 0)
	{
		factor *= get_value("rCC");

	} // Child (2) - Mother (0)
	else if (get_infer("rMC") && ((H->h_residents[idx].i_age == 0 && H->h_residents[i].i_age == 1 && H->h_residents[i].i_sex == 1) || (H->h_residents[i].i_age == 0 && H->h_residents[idx].i_age == 1 && H->h_residents[idx].i_sex == 1)))
	{
		factor *= get_value("rMC");

	} // Child (2) - Father (1)
	else if (get_infer("rFC") && ((H->h_residents[idx].i_age == 0 && H->h_residents[i].i_age == 1 && H->h_residents[i].i_sex == 2) || (H->h_residents[i].i_age == 0 && H->h_residents[idx].i_age == 1 && H->h_residents[idx].i_sex == 2)))
	{
		factor *= get_value("rFC");
	}

	return factor;
}

/* Model the susceptibility within household */
double Model::susceptibility(Household *H, size_t idx)
{
	double factor = 1.0;

	// Suceptibility of different vaccination classes
	if (get_infer("rU") && H->h_residents[idx].i_vacc_status == 1)
	{
		factor *= get_value("rU");
	}

	if (get_infer("r1V") && H->h_residents[idx].i_vacc_status == 2 && H->h_residents[idx].i_vacc_delay < get_value("vacc_onset"))
	{
		factor *= get_value("r1V");
	}

	if (get_infer("r1V90") && H->h_residents[idx].i_vacc_status == 2 && H->h_residents[idx].i_vacc_delay >= get_value("vacc_onset"))
	{
		factor *= get_value("r1V90");
	}

	if (get_infer("r2V") && H->h_residents[idx].i_vacc_status == 3 && H->h_residents[idx].i_vacc_delay < get_value("vacc_onset"))
	{
		factor *= get_value("r2V");
	}

	// Suceptibility of children
	if (get_infer("rC") && H->h_residents[idx].i_age == 0)
	{
		factor *= get_value("rC");
	}

	// Suceptibility of isolated people
	if (get_infer("rQ") && H->h_residents[idx].i_isolation)
	{
		factor *= get_value("rQ");
	}

	// Suceptibility of recovered people
	if (get_infer("rR") && H->h_residents[idx].i_recovered && (H->h_residents[idx].i_vacc_status == 1 || H->h_residents[idx].i_vacc_status == 2))
	{
		factor *= get_value("rR");
	}

	return factor;
}

/* Model the transmissability within household */
double Model::transmissability(Household *H)
{
	double factor = 1.0;

	// Instantaneous per capita transmission rate
	factor *= get_value("beta");

	// If transmission hazard depends on household size
	if (get_infer("delta"))
	{
		factor /= pow(H->h_size / get_value("typical_size"), get_value("delta"));
	}

	return factor;
}

/* Poisson point process
void Model::simulate_data(timetype delta_t, std::mt19937_64 &rng)
{
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	timetype infec_time;
	size_t res;

	// initialize all residents to healthy
	H->h_healthy_list.resize(0);
	H->h_symptomatic_list.resize(0);
	for (size_t res = 0; res < H->h_size; res++)
	{
		set_status(res, HEALTHY);
		h_residents[res].init_augm("incub_period", rng);
		set_augm("onset_time", res, INFINITY);
		h_healthy_list.push_back(res);
	}
	// initialize h_infec_time
	h_infec_time = INFINITY;
	h_first_case = h_size;

	// FIXME : upper limit
	for (timetype infec_time = 0; infec_time < 400; infec_time += delta_t)
	{
		for (size_t i = 0; i < h_healthy_list.size(); i++)
		{
			res = h_healthy_list[i];
			set_augm("onset_time", res, infec_time + get_augm("incub_period", res));

			// FIXME : the survival period inside here is wrong
			// it should be just the interval
			if (unif(rng) < std::exp(log_infec_proba(res)))
			{
				set_status(res, SYMPTOMATIC);
				h_symptomatic_list.push_back(res);
				h_healthy_list.erase(h_healthy_list.begin() + i);
				// FIXME : h_first and h_infec to be updated
			}
			else
			{
				set_augm("onset_time", res, INFINITY);
			}
		}
	}

	// outside
	// renormalize to have the same infection period of dataset
}
*/
