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
	int status,
	int indid,
	int age,
	int sex,
	bool isolation,
	bool recovered,
	int vacc_status,
	timetype study_period,
	timetype vacc_delay,
	timetype onset_time)
{
	// Status
	i_status = status;
	// Data variables
	i_indid = indid;
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
						  double value_0,
						  kwargs value_range)
{
	Augmented this_augm(new_value, model_name, model_args_ext, strategy_name, strategy_args, value_0, value_range);
	i_augm.emplace(augmID, this_augm);
	i_augm_list.push_back(augmID);

	if (model_name != UNDEF)
	{
		i_augm_with_model.push_back(augmID);
	}
	if (strategy_name != UNDEF)
	{
		i_augm_to_be_inferred.push_back(augmID);
	}
}

double Individual::augm_model_log_pdf(std::string augmID)
{
	double new_value = i_augm.at(augmID).get_value();

	return i_augm.at(augmID).model_log_pdf(new_value);
}

double Individual::augm_model_log_ratio(std::string augmID, double old_value, double new_value)
{
	return i_augm.at(augmID).model_log_ratio(old_value, new_value);
}

double Individual::augm_propose(std::string augmID, std::mt19937_64 &rng)
{
	return i_augm.at(augmID).propose(rng);
}

double Individual::augm_propose_log_ratio(std::string augmID, double old_value, double new_value)
{
	return i_augm.at(augmID).propose_log_ratio(old_value, new_value);
}

/////////////////////////////////////////////////////
//  Household : Class
/////////////////////////////////////////////////////

// UPDATE : circular / temporary

// Add individual to household
void Household::add_individual(Individual this_indiv)
{
	h_residents.push_back(this_indiv);
	set_status(h_size, h_residents[h_size].get_status());
	h_size++;
}

void Household::set_status(size_t idx, int status)
{
	h_residents[idx].set_status(status);
	if (status == HEALTHY)
	{
		// Update the list of healthy residents within the household
		h_healthy_list.push_back(idx);
	}
	else
	{
		// Update the list of symptomatic residents within the household
		h_symptomatic_list.push_back(idx);
	}
}