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
	timetype study_period,
	timetype vacc_delay,
	timetype onset_time)
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

void Individual::add_augm(std::string augmID,
						  double new_value,
						  std::string model_name,
						  kwargs_ext model_args_ext,
						  std::string strategy_name,
						  kwargs strategy_args,
						  kwargs value_range)
{
	Augmented this_augm(new_value, model_name, model_args_ext, strategy_name, strategy_args, value_range);
	i_augm.emplace(augmID, this_augm);
	i_augm_list.push_back(augmID);

	if (model_name != UNDEF)
	{
		i_augm_with_model_list.push_back(augmID);
	}
}

double Individual::augm_propose(std::string augmID, std::mt19937_64 &rng)
{
	// NOTE : I choose to propose augmented values around the initial value
	double old_value = i_augm.at(augmID).get_value_0();
	// FIXME: this should be a user variable, probably in metropolis in order to propose from 0 or not

	return i_augm.at(augmID)._base_.propose(rng, old_value);
}

double Individual::augm_model_log_pdf(std::string augmID)
{
	double new_value = i_augm.at(augmID).get_value();

	return i_augm.at(augmID).model_log_pdf(new_value);
}

double Individual::augm_model_log_ratio(std::string augmID, double new_value)
{
	// FIXME : I choose to propose augmented values around the initial value
	double old_value = i_augm.at(augmID).get_value_0();

	// FIXME : what if there's no model? i.e. uniform always 0?

	return i_augm.at(augmID).model_log_ratio(old_value, new_value);
}

double Individual::augm_propose_log_ratio(std::string augmID)
{
	// FIXME : I choose to propose augmented values around the initial value

	return i_augm.at(augmID)._base_.propose_log_ratio();
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

// Initialize augmented variables and infectivity profiles
void Household::init_augm_data(std::mt19937_64 &rng)
{
	// loop over residents within the house
	for (size_t i = 0; i < h_size; i++)
	{
		h_residents[i].init_augm_data(rng);
	}
	init_household_status();
}

// Update list of healthy and symptomatic residents
void Household::init_household_status()
{
	h_healthy_list.resize(0);
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
		}
	}
}

// Update household first case
void Household::update_household_first_case()
{
	timetype t_infec, new_t0 = INFINITY;
	size_t new_case0 = h_size; // NOTE: init to impossible values

	for (size_t i : h_symptomatic_list)
	{
		// Find smaller infection time
		t_infec = h_residents[i].get_infec_time();
		new_t0 = std::min(t_infec, new_t0);

		// Identify first case
		if (t_infec == new_t0)
		{
			new_case0 = i;
		}
	}

	h_infec_time = new_t0;
	h_first_case = new_case0;
}

// ----------------------------
//  Methods : Aliases
// ----------------------------

// Change not healthy individual augmented variables
void Household::set_augm(std::string augmID, size_t idx, double new_value)
{
	h_residents[idx].set_augm(augmID, new_value);
	update_household_first_case();	  // Update information of first cases
	update_infectivity_profiles(idx); // Update infectivity profiles realted to idx
}

// ----------------------------
//  Methods : Computations
// ----------------------------

/* Initialize infectivity profiles */
void Household::init_infectivity_profiles(Onset_to_Transm *infectivity_profile)
{
	ptr_to_infec_profile = infectivity_profile;

	// initialize 2-dim vectors h_cumul_profile and h_inst_profile
	h_cumul_profile.resize(h_size);
	h_inst_profile.resize(h_size);

	// loop over all residents
	for (size_t i = 0; i < h_size; i++)
	{
		h_residents[i].ptr_to_infec_profile = infectivity_profile;
		h_cumul_profile[i].resize(h_size);
		h_inst_profile[i].resize(h_size);
	}
}

/* Update infectivity profiles given a NOT healthy resident idx */
void Household::update_infectivity_profiles(size_t idx)
{
	timetype t_infec_idx = get_infec_time(idx), t_infec_res;

	// Update cumulative and instantaneous person-to-person transmission rates
	// NOTE : h_cumul_profile[idx][res] : from `idx` to `res`

	// loop over all non healthy residents
	for (size_t res : h_symptomatic_list)
	{
		t_infec_res = get_infec_time(res);
		if (t_infec_idx < t_infec_res)
		{ // `idx` : infector, `res` : infectee
			h_cumul_profile[idx][res] = h_residents[idx].infectivity_cumul(t_infec_res);
			h_inst_profile[idx][res] = h_residents[idx].infectivity_inst(t_infec_res);
		}
		else if (t_infec_res < t_infec_idx)
		{ // `res` : infector, `idx` : infectee
			h_cumul_profile[res][idx] = h_residents[res].infectivity_cumul(t_infec_idx);
			h_inst_profile[res][idx] = h_residents[res].infectivity_inst(t_infec_idx);
		}
		// NOTE : else is the case idx == res
	}

	// loop over all healthy residents
	// UPDATE : if study_period is equal for all residents, than this is constant for all
	for (size_t res : h_healthy_list)
	{
		t_infec_res = get_infec_time(res); // == study_period
		h_cumul_profile[idx][res] = h_residents[idx].infectivity_cumul(t_infec_res);
	}
}
