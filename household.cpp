// household.cpp //

#include "household.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Individual : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

Individual::Individual(
	int indid,
	int status,
	int age,
	int sex,
	bool isolation,
	bool recovered,
	int vacc_status,
	double study_period,
	double vacc_delay,
	double onset_time)
{
	// Data variables
	i_indid = indid;
	i_status = status;
	i_age = age;
	i_sex = sex;
	i_isolation = isolation;
	i_recovered = recovered;
	i_vacc_status = vacc_status;
	i_study_period = study_period;
	i_vacc_delay = vacc_delay;
	i_onset_time = onset_time;
}

double Individual::propose_augm(std::string augmID, std::mt19937_64 &rand_gen)
{
	// NOTE : I choose to propose augmented values around the initial value
	double old_value = i_augm.at(augmID).get_value_0();
	// FIXME: this should be a user variable, probably in metropolis

	return i_augm.at(augmID).propose(rand_gen, old_value);
}

double Individual::log_ratio_prior_augm(std::string augmID, double new_value)
{
	// NOTE : I choose to propose augmented values around the initial value
	double old_value = i_augm.at(augmID).get_value_0();
	
	return i_augm.at(augmID).prior_log_ratio(old_value, new_value);
}

double Individual::log_ratio_proposal_augm(std::string augmID, double new_value)
{
	// NOTE : I choose to propose augmented values around the initial value
	double old_value = i_augm.at(augmID).get_value_0();
	
	return i_augm.at(augmID).propose_log_ratio(old_value, new_value);
}

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
// UPDATE: how to formalize the model dependency on augmented variable?

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
		if (get_status(i) == HEALTHY)
		{
			// Update the list of healthy residents within the household
			h_healthy_list.push_back(i);
		}
		else
		{
			// Update the list of symptomatic residents within the household
			h_symptomatic_list.push_back(i);

			// Find smaller infection time
			new_t0 = std::min(h_residents[i].get_infec_time(), new_t0);

			// Identify first case
			if (h_residents[i].get_infec_time() == new_t0)
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
void Household::set_augm(std::string augmID, size_t idx, double new_incub_period)
{
	h_residents[idx].set_augm(augmID, new_incub_period);
	update_infectivity_profiles(idx);
}

// Propose new value for symptom onset time
double Household::propose_augm(std::string augmID, size_t idx, std::mt19937_64 &rand_gen)
{
	return h_residents[idx].propose_augm(augmID, rand_gen);
}

// Check if the augmented variable is within range
bool Household::valid_augm(std::string augmID, size_t idx, double new_value)
{
	return h_residents[idx].valid_augm(augmID, new_value);
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
			t_transm = get_infec_time(res);

			if (get_infec_time(idx) < t_transm)
			{ // NOTE : this ensures idx != res

				// Cumulative transmission rate for all household contacts
				h_cumul_rate[idx][res] = h_residents[idx].infectivity_cumul(t_transm);

				// Instantaneous transmission rate for secondary cases
				if (get_status(res) != HEALTHY)
				{ // resident `res` got infected
					h_inst_rate[idx][res] = h_residents[idx].infectivity_inst(t_transm);
				}
			}
		}
	}

	// --------------------------------------------------------
	// resident `idx` is considered as a potential infectee
	// --------------------------------------------------------

	study_period = h_residents[idx].i_study_period;

	// time of transmission to `idx`
	t_transm = get_infec_time(idx);

	// loop over all possible infectors, i. e. resident `j`
	for (size_t res : h_symptomatic_list)
	{
		if (get_infec_time(res) < t_transm)
		{ // NOTE : this ensures idx != res

			// Cumulative transmission rate for all household contacts
			h_cumul_rate[res][idx] = h_residents[res].infectivity_cumul(t_transm);

			// Instantaneous transmission rate for secondary cases
			if (get_status(idx) != HEALTHY)
			{ // resident `idx` got infected
				h_inst_rate[res][idx] = h_residents[res].infectivity_inst(t_transm);
			}
		}
	}
}
