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
void Individual::init_augm_data(std::mt19937_64& rand_gen)
{
	bool is_sick = (i_status != HEALTHY);
	double value;

	// FIXME : args change with mcmc : I need pointers, not values!

	/* onset_time */
	{
		std::string augmID = "onset_time";
		double init_value = (is_sick == true) ? i_onset_time : i_study_period;
		std::string strategy = "normal";
		kwargs args = {
			{"mean", 0.0},
			{"sd", glob_params.get_value("onset_err")}};
		kwargs range = {{"low", -MAX_INCUB_PERIOD}, {"up", i_study_period}};

		Augmented_Variable this_augm(init_value, is_sick, strategy, args, range);
		i_augm.emplace(augmID, this_augm);
		i_augm_list.push_back(augmID);
	}

	/* incub_period */
	{
		std::string augmID = "incub_period";
		double init_value = 0.;
		std::string strategy = "gamma";
		kwargs args = {
			{"shape", glob_params.get_value("incub_shape")},
			{"scale", glob_params.get_value("incub_scale")}};
		kwargs range = {{"low", 0.0}, {"up", MAX_INCUB_PERIOD}};

		Augmented_Variable this_augm(init_value, is_sick, strategy, args, range);
		i_augm.emplace(augmID, this_augm);
		i_augm_list.push_back(augmID);
	}

	// initialize for unhealthy
	if (is_sick == true)
	{
		for (std::string augmID : i_augm_list)
		{
			value = propose_augm(augmID, rand_gen);
			while (valid_augm(augmID, value) == false)
			{
				value = propose_augm(augmID, rand_gen);
			}
			set_augm(augmID, value);
		}
	}
}

/* Time of infection */
double Individual::get_infec_time() const
{
	return get_augm("onset_time") - get_augm("incub_period");
}

// ----------------------------
//  Methods : Infectivity Profiles
// ----------------------------

// FIXME : these are probabilities, but we want rates!!!

double Individual::infectivity_inst(double t_transm)
{
	double output = 0.0, tost;
	kwargs args = {
		{"alpha", glob_params.get_value("tost_alpha")},
		{"mu", glob_params.get_value("tost_mu")},
		{"sigma", glob_params.get_value("tost_sigma")},
		{"tau", glob_params.get_value("tost_tau")}};

	tost = t_transm - get_augm("onset_time");
	output = t_onset_sympt_to_transm_pdf(tost, get_augm("incub_period"), args);

	return output;
}

double Individual::infectivity_cumul(double t_transm)
{
	double output = 0.0, tost;
	kwargs args = {
		{"alpha", glob_params.get_value("tost_alpha")},
		{"mu", glob_params.get_value("tost_mu")},
		{"sigma", glob_params.get_value("tost_sigma")},
		{"tau", glob_params.get_value("tost_tau")}};

	tost = t_transm - get_augm("onset_time");
	output = t_onset_sympt_to_transm_cdf(tost, get_augm("incub_period"), args);

	return output;
}

/////////////////////////////////////////////////////
//  Household : Methods (likelihood)
/////////////////////////////////////////////////////

/* logarithm of the incubation probability */
double Household::log_incub_proba(size_t idx)
{
	double x = get_augm("incub_period", idx);
	kwargs args = {
		{"shape", glob_params.get_value("incub_shape")},
		{"scale", glob_params.get_value("incub_scale")}};

	return log_gamma_pdf(x, args);
	// FIXME : this should come from resident augm prior pdf

	// return h_residents[idx].log_incub_period_pdf(x);
}

// FIXME : these two below can be compacted to avoid errors

double Household::model_risk_beta(const std::vector<std::vector<double>>& rate_vector, size_t idx)
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

/* logarithm of the infection probability */
double Household::log_infec_proba(size_t idx)
{
	double beta = 0.0, output = 0.0;
	double t_transm = get_infec_time(idx);

	beta = model_risk_beta(h_inst_rate, idx);
	output = model_susceptibility(idx) * (model_community_risk_inst(t_transm) + beta);
	output = std::log(1. - std::exp(-output));

	output = output + log_surv_proba(idx);

	return output;
}

/* Survival probability from time t0 of resident idx */
double Household::log_surv_proba(size_t idx)
{
	double beta = 0.0, output = 0.0;
	double t0 = h_infec_time;
	double t_transm = get_infec_time(idx);

	beta = model_risk_beta(h_cumul_rate, idx);
	output = model_susceptibility(idx) * (model_community_risk_cumul(t0, t_transm) + beta);
	output *= (-1.);

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
		// FIXME : contribution of augmented variables priors
		loglike_house += log_incub_proba(res);

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
double Household::model_community_risk_inst(double t1)
{
	double output = 0.0;
	output = lin_idx_interp(t1, Inst_Community_Infec);
	output *= glob_params.get_value("alpha");

	return output;
}

// Cumulative risk of infection from community
double Household::model_community_risk_cumul(double t0, double t1)
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