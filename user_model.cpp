// model.cpp //

#include "household.hpp"

using namespace std;

// Shared Parameter Dictionary
extern std::vector<double> Cumul_Community_Infec;
extern std::vector<double> Inst_Community_Infec;

/////////////////////////////////////////////////////
//  Individual : Methods (likelihood)
/////////////////////////////////////////////////////

// ----------------------------
//  Methods : Augmented Variables
// ----------------------------

// List of augmented variables
std::vector<std::string> List_Augm_Var = {"onset_time", "incub_period"};

/* Initialize augmented variables */
// UPDATE: this needs to become more user friendly
void Individual::init_augm_data(std::mt19937_64 &rng)
{
	// --- onset_time --- //
	{
		std::string augmID = "onset_time";
		double init_value = (i_status == HEALTHY) ? i_study_period : i_onset_time;

		/*
		std::string model_name = UNDEF;
		kwargs_ext model_args = {};
		*/
		// FIXME : how to deal with these cases!?
		std::string model_name = "unif";
		kwargs_ext model_args = {
			{"low", glob_params.get_address("onset_low")},
			{"up", glob_params.get_address("onset_up")}};

		std::string metropolis_strategy = (i_status == HEALTHY) ? UNDEF : "normal";

		kwargs metropolis_args = {
			{"mean", 0.0},
			{"sd", glob_params.get_value("onset_err")}};

		kwargs range = {{"low", glob_params.get_value("onset_low")}, {"up", glob_params.get_value("onset_up")}};

		add_augm(augmID, init_value, model_name, model_args, metropolis_strategy, metropolis_args, range);
	}

	// --- incub_period --- //
	{
		std::string augmID = "incub_period";
		double init_value = 0.;

		std::string model_name = (i_status == HEALTHY) ? UNDEF : "gamma";
		kwargs_ext model_args = {
			{"shape", glob_params.get_address("incub_shape")},
			{"scale", glob_params.get_address("incub_scale")}};

		// FIXME : waiting for prior proposal
		std::string metropolis_strategy = (i_status == HEALTHY) ? UNDEF : "normal";
		kwargs metropolis_args = {
			{"mean", glob_params.get_value("incub_shape")},
			{"sd", glob_params.get_value("incub_scale")}};

		kwargs range = {{"up", MAX_INCUB_PERIOD}};

		add_augm(augmID, init_value, model_name, model_args, metropolis_strategy, metropolis_args, range);
		if (model_name != UNDEF)
		{
			set_augm(augmID, augm_model_rand(augmID, rng));
		}
	}
}

/* Time of infection */
timetype Individual::get_infec_time() const
{
	return get_augm("onset_time") - get_augm("incub_period");
}

// ----------------------------
//  Methods : Infectivity Profiles
// ----------------------------

double Individual::infectivity_inst(timetype t_transm)
{
	timetype tost = t_transm - get_augm("onset_time");
	return ptr_to_infec_profile->pdf(tost, get_augm("incub_period"));
}

double Individual::infectivity_cumul(timetype t_transm)
{
	timetype tost = t_transm - get_augm("onset_time");
	return ptr_to_infec_profile->cdf(tost, get_augm("incub_period"));
}

/////////////////////////////////////////////////////
//  Household : Methods (likelihood)
/////////////////////////////////////////////////////

double Household::model_risk_beta(const std::vector<std::vector<double>> &rate_vector, size_t idx)
{
	double beta_i = 0.0, beta = 0.0;

	// Risk of infection from other household members
	for (size_t res : h_symptomatic_list)
	{
		if (idx != res)
		{
			beta_i = rate_vector[res][idx];
			beta_i *= model_contact_modulation(res, idx);
			beta += beta_i;
		}
	}

	beta *= model_transmissability();

	return beta;
}

/* Survival probability from time t0 of resident idx */
double Household::log_surv_proba(size_t idx)
{
	double beta = 0.0, output = 0.0;
	timetype t0 = h_infec_time;
	timetype t_transm = get_infec_time(idx);

	beta = model_risk_beta(h_cumul_profile, idx);
	output = model_susceptibility(idx) * (model_community_risk_cumul(t0, t_transm) + beta);
	output *= (-1.);

	return output;
}

/* logarithm of the infection probability */
double Household::log_infec_proba(size_t idx)
{
	double beta = 0.0, output = 0.0;
	timetype t_transm = get_infec_time(idx);

	beta = model_risk_beta(h_inst_profile, idx);
	output = model_susceptibility(idx) * (model_community_risk_inst(t_transm) + beta);
	output = std::log(1. - std::exp(-output));

	output = output + log_surv_proba(idx);

	return output;
}

/* Compute the log-likelihood within household */
double Household::compute_loglike()
{
	double loglike_house = 0.0;

	// Loop over all healthy residents // NOTE: no augmented variables

	for (size_t res : h_healthy_list)
	{
		loglike_house += log_surv_proba(res);
	}
	// Loop over all symptomatic residents
	for (size_t res : h_symptomatic_list)
	{
		// Contribution of model associated to augmented variables
		for (std::string augmID : h_residents[res].i_augm_with_model_list)
		{
			loglike_house += augm_model_log_pdf(augmID, res);
		}

		// Contribution of secondary cases
		if (res != h_first_case)
		{
			loglike_house += log_infec_proba(res);
		}
	}

	return loglike_house;
}

/////////////////////////////////////////////////////
//  Household : Methods (modulators)
/////////////////////////////////////////////////////

// UPDATE : this could be interpolated instead

// Instaneous risk of infection from community
double Household::model_community_risk_inst(timetype t1)
{
	double output = 0.0;
	output = lin_idx_interp(t1, Inst_Community_Infec);
	output *= glob_params.get_value("alpha");

	return output;
}

// Cumulative risk of infection from community
double Household::model_community_risk_cumul(timetype t0, timetype t1)
{
	double output = 0.0;
	output = lin_idx_interp(t1, Cumul_Community_Infec) - lin_idx_interp(t0, Cumul_Community_Infec);
	output *= glob_params.get_value("alpha");

	return output;
}

// Relative infectivity of symptomatic child cases versus symptomatic adult cases
double Household::model_contact_modulation(size_t infectee, size_t idx)
{
	double factor = 1.0;
	size_t i = infectee;

	// NOTE : this call is missleading : use only if to be inferred...

	// Child (2) - Child (2)
	if (glob_params.get_infer("rCC") == true && h_residents[idx].i_age == 0 && h_residents[i].i_age == 0)
	{
		factor *= glob_params.get_value("rCC");

	} // Child (2) - Mother (0)
	else if (glob_params.get_infer("rMC") == true && ((h_residents[idx].i_age == 0 && h_residents[i].i_age == 1 && h_residents[i].i_sex == 1) || (h_residents[i].i_age == 0 && h_residents[idx].i_age == 1 && h_residents[idx].i_sex == 1)))
	{
		factor *= glob_params.get_value("rMC");

	} // Child (2) - Father (1)
	else if (glob_params.get_infer("rFC") == true && ((h_residents[idx].i_age == 0 && h_residents[i].i_age == 1 && h_residents[i].i_sex == 2) || (h_residents[i].i_age == 0 && h_residents[idx].i_age == 1 && h_residents[idx].i_sex == 2)))
	{
		factor *= glob_params.get_value("rFC");
	}

	return factor;
}

double Household::model_susceptibility(size_t idx)
{
	double factor = 1.0;

	// Suceptibility of different vaccination classes
	if (glob_params.get_infer("rU") == true && h_residents[idx].i_vacc_status == 1)
	{
		factor *= glob_params.get_value("rU");
	}

	if (glob_params.get_infer("r1V") == true && h_residents[idx].i_vacc_status == 2 && h_residents[idx].i_vacc_delay < VACC_ONSET)
	{
		factor *= glob_params.get_value("r1V");
	}

	if (glob_params.get_infer("r1V90") == true && h_residents[idx].i_vacc_status == 2 && h_residents[idx].i_vacc_delay >= VACC_ONSET)
	{
		factor *= glob_params.get_value("r1V90");
	}

	if (glob_params.get_infer("r2V") == true && h_residents[idx].i_vacc_status == 3 && h_residents[idx].i_vacc_delay < VACC_ONSET)
	{
		factor *= glob_params.get_value("r2V");
	}

	if (glob_params.get_infer("r3V") == true && h_residents[idx].i_vacc_status == 4)
	{
		factor *= glob_params.get_value("r3V");
	}

	// Suceptibility of children
	if (glob_params.get_infer("rC") == true && h_residents[idx].i_age == 0)
	{
		factor *= glob_params.get_value("rC");
	}

	// Suceptibility of isolated people
	if (glob_params.get_infer("rQ") == true && h_residents[idx].i_isolation == true)
	{
		factor *= glob_params.get_value("rQ");
	}

	// Suceptibility of recovered people
	if (glob_params.get_infer("rR") == true && h_residents[idx].i_recovered == true && (h_residents[idx].i_vacc_status == 1 || h_residents[idx].i_vacc_status == 2))
	{
		factor *= glob_params.get_value("rR");
	}

	return factor;
}

double Household::model_transmissability()
{
	double factor = 1.0;

	// Instantaneous per capita transmission rate
	factor *= glob_params.get_value("beta");

	// If transmission hazard depends on household size
	if (glob_params.get_infer("delta") == true)
	{
		factor /= pow(h_size / glob_params.get_value("typical_size"), glob_params.get_value("delta"));
	}

	return factor;
}