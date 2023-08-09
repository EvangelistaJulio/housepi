// household.cpp //

#include "household.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Household : Class
/////////////////////////////////////////////////////

// Add individual to household
void Household::add_individual(Individual this_indiv)
{
	h_residents.push_back(this_indiv);
	h_size++;
}

// logarithm of the pdf of the incubation period
double Household::log_incub_period_pdf(size_t idx)
{
	double x = get_augm_incub_period(idx);
	return h_residents[idx].log_incub_period_pdf(x);
}

// Initialize augmented variables and infectivity profiles
void Household::init_augm_data(std::mt19937_64 &rand_gen)
{
	// loop over residents within the house
	for (size_t i = 0; i < h_size; i++)
	{
		h_residents[i].init_augm_data(rand_gen);
	}
	update_household_status();
	init_infectivity_profiles();
}

// Update list of symptomatic residents and first case
// NOTE: repeating these operation at fixed status is expensive
void Household::update_household_status()
{
	double new_t0 = INFINITY;
	size_t new_case0 = h_size; // NOTE: init to impossible value

	h_symptomatic_list.resize(0);
	for (size_t i = 0; i < h_size; i++)
	{
		if (get_status(i) != HEALTHY)
		{
			// Update the list of symptomatic residents within the household
			h_symptomatic_list.push_back(i);

			// Find smaller infection time
			new_t0 = std::min(h_residents[i].get_augm_infec_time(), new_t0);

			// Identify first case
			if (h_residents[i].get_augm_infec_time() == new_t0)
			{
				new_case0 = i;
			}
		}
	}

	h_infec_time = new_t0;
	h_first_case = new_case0;
}

// ----------------------------
//  Methods : Aliases
// ----------------------------

// Change individual infection time
void Household::set_augm_incub_period(size_t idx, double new_incub_period)
{
	h_residents[idx].set_augm_incub_period(new_incub_period);
	update_infectivity_profiles(idx);
}

// Change symptom onset
void Household::set_augm_onset_time(size_t idx, double new_onset_time)
{
	h_residents[idx].set_augm_onset_time(new_onset_time);
	update_infectivity_profiles(idx);
}

// Propose new value for symptom onset time
double Household::propose_augm_onset_time(size_t idx, std::mt19937_64 &rand_gen)
{
	return h_residents[idx].propose_augm_onset_time(rand_gen);
}

// Propose new value for infection time
double Household::propose_augm_incub_period(size_t idx, std::mt19937_64 &rand_gen)
{
	return h_residents[idx].propose_augm_incub_period(rand_gen);
}

// ----------------------------
//  Methods : Computations
// ----------------------------

/* Initialize infectivity profiles */
void Household::init_infectivity_profiles()
{
	// initialize 2-dim vectors h_cumul_rate and h_inst_rate
	h_cumul_rate.resize(h_size);
	h_inst_rate.resize(h_size);

	// loop over all residents
	for (size_t i = 0; i < h_size; i++)
	{
		h_cumul_rate[i].resize(h_size);
		h_inst_rate[i].resize(h_size);
	}

	// Compute cumulative and instantaneous person-to-person transmission rates
	for (size_t i : h_symptomatic_list)
	{
		update_infectivity_profiles(i);
	}
}

/* Update infectivity profiles */
void Household::update_infectivity_profiles(size_t idx)
{
	double t_transm, study_period;

	// Update information of symptomatic and first cases
	update_household_status();

	// Update cumulative and instantaneous person-to-person transmission rates
	// NOTE : h_cumul_rate[idx][res] : from `idx` to `res`

	// --------------------------------------------------------
	// resident `idx` is considered as a potential infector
	// --------------------------------------------------------
	if (get_status(idx) != HEALTHY)
	{
		// loop over all possible infecteee, i. e. resident `res`
		for (size_t res = 0; res < h_size; res++)
		{
			study_period = h_residents[res].i_study_period;

			// time of transmission to resident `r`
			t_transm = get_augm_infec_time(res);

			if (get_augm_infec_time(idx) < t_transm)
			{ // NOTE : this ensures idx != res

				// Cumulative transmission rate for all household contacts
				h_cumul_rate[idx][res] = h_residents[idx].infec_cumul(t_transm);

				// Instantaneous transmission rate for secondary cases
				if (get_status(res) != HEALTHY)
				{ // resident `res` got infected
					h_inst_rate[idx][res] = h_residents[idx].infec_inst(t_transm);
				}
			}
		}
	}

	// --------------------------------------------------------
	// resident `idx` is considered as a potential infectee
	// --------------------------------------------------------

	study_period = h_residents[idx].i_study_period;

	// time of transmission to `idx`
	t_transm = get_augm_infec_time(idx);

	// loop over all possible infectors, i. e. resident `j`
	for (size_t res : h_symptomatic_list)
	{
		if (get_augm_infec_time(res) < t_transm)
		{ // NOTE : this ensures idx != res

			// Cumulative transmission rate for all household contacts
			h_cumul_rate[res][idx] = h_residents[res].infec_cumul(t_transm);

			// Instantaneous transmission rate for secondary cases
			if (get_status(idx) != HEALTHY)
			{ // resident `idx` got infected
				h_inst_rate[res][idx] = h_residents[res].infec_inst(t_transm);
			}
		}
	}
}

/* logarithm of the infection probability */
double Household::log_infec_proba(size_t idx)
{
	double beta_i = 0.0, beta = 0.0, output = 0.0;
	double t_transm = get_augm_infec_time(idx);

	// Instantaneous risk of infection from other household members
	for (size_t res : h_symptomatic_list)
	{
		if (idx != res)
		{
			beta_i = h_inst_rate[res][idx];
			beta_i *= model_contact_modulation(res, idx);
			beta += beta_i;
		}
	}

	beta *= model_transmissability();
	output = model_susceptibility(idx) * (model_community_risk_inst(t_transm) + beta);
	output = std::log(1. - std::exp(-output));

	return output;
}

/* Survival probability from time t0 of resident idx */
double Household::log_surv_proba(size_t idx, double t0)
{
	double output = 0.0, beta = 0.0, beta_i = 0.0;
	double t_transm = get_augm_infec_time(idx);

	// Cumulative risk of infection from other household residents
	for (size_t res : h_symptomatic_list)
	{
		if (idx != res)
		{
			beta_i = h_cumul_rate[res][idx];
			beta_i *= model_contact_modulation(res, idx);
			beta += beta_i;
		}
	}

	beta *= model_transmissability();
	output = model_susceptibility(idx) * (model_community_risk_cumul(t0, t_transm) + beta);
	output *= (-1.);

	return output;
}

/* Compute the log-likelihood within household */
double Household::compute_loglike()
{
	double loglike_house = 0.0, loglike_res = 0.0;

	// Loop over all residents
	for (size_t res = 0; res < h_size; res++)
	{
		if (h_residents[res].i_status == HEALTHY)
		{ // Contribution of healthy cases
			loglike_res = log_surv_proba(res, h_infec_time);
		}
		else if (res == h_first_case)
		{ // Contribution of the first case
			loglike_res = log_incub_period_pdf(res);
		}
		else
		{ // Contribution of secondary cases
			loglike_res = log_incub_period_pdf(res);
			loglike_res += log_surv_proba(res, h_infec_time);
			loglike_res += log_infec_proba(res);
		}
		loglike_house += loglike_res;
	}

	return loglike_house;
}