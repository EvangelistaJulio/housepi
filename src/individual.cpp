// individual.cpp //

#include "individual.h"

#define MAX_NUMBER_ATTEMPTS 10

////////////////////////////////////////////////////////////////////////
// /////////////////// Individual : Constructor ///////////////////// //
////////////////////////////////////////////////////////////////////////

Individual::Individual(sdargs new_data)
{
	index = new_data.at("index");
	status = new_data.at("status");
	data = new_data;
}

/* get_data : bottleneck if data are called every iteration. Deprecated outisde of pre-processing.*/
double Individual::get_data(std::string dataID) const
{
	if (data.find(dataID) == data.end())
	{
		std::string err_text = "Data column `" + dataID + "` unknown.";
		exit_with_failure_error(err_text);
	}
	return data.at(dataID);
}

//=====Rate Methods (e.g. Infectivity, Susceptibility, ...)=====//

/* Assign a model rate */
void Individual::assign_rate_per_group(std::string rateID, std::string parID)
{
	if (_which_rate.find(rateID) == _which_rate.end())
	{
		// if rateID was undeclared, add it
		_which_rate.emplace(rateID, parID);
	}
	else
	{
		// rateID was previously declared
		std::string err_text = "Found individual with multiple definitions for group rate `" + rateID + "`.";
		exit_with_failure_error(err_text);
	}
}

/* Link the rate value with a pointer */
void Individual::_set_rate(std::string rateID, const double *x)
{
	_rate_ptr.emplace(rateID, x);
}

/* Return value of model rate */
double Individual::model_rate(std::string rateID) const
{
	return *_rate_ptr.at(rateID);
}

//=====Augmented Variables=====//

/* Init augmented variable within a random initial value */
void Individual::init_augm(std::string augmID, std::mt19937_64 &rng)
{
	int num_attempts = 0;
	double init_value = augm_model_rand(augmID, rng);

	while ((!augm_valid(augmID, init_value)))
	{
		if (num_attempts >= MAX_NUMBER_ATTEMPTS)
		{
			std::string err_text = "Number of attempts exceed in random initialization of `" + augmID + "`.";
			exit_with_failure_error(err_text);
		}
		init_value = augm_model_rand(augmID, rng);
		num_attempts++;
	}
	set_augm(augmID, init_value);
}

/* Add an augmented variable to the individual. */
void Individual::add_augm(std::string augmID, Augmented new_augm)
{
	if (augm_var.find(augmID) == augm_var.end())
	{
		augm_var.emplace(augmID, new_augm);
	}
	else
	{
		std::string err_text = "Augmented variable `" + augmID + "` has multiple definitions.";
		exit_with_failure_error(err_text);
	}
}

/* Add an augmented variable to the individual, plus initialize it. */
void Individual::add_augm(std::string augmID, Augmented new_augm, std::mt19937_64 &rng)
{
	add_augm(augmID, new_augm);
	init_augm(augmID, rng);
}

double Individual::augm_model_log_pdf(std::string augmID)
{
	double new_val = augm_var.at(augmID).get_value();
	return augm_var.at(augmID).model_log_pdf(new_val);
}

double Individual::augm_model_log_ratio(std::string augmID, double old_val, double new_val)
{
	return augm_var.at(augmID).model_log_ratio(old_val, new_val);
}

double Individual::augm_propose(std::string augmID, std::mt19937_64 &rng)
{
	return augm_var.at(augmID).propose(rng);
}

double Individual::augm_propose_log_ratio(std::string augmID, double old_val, double new_val)
{
	return augm_var.at(augmID).propose_log_ratio(old_val, new_val);
}

////////////////////////////////////////////////////////////////////////
// /////////////////////////// Household //////////////////////////// //
////////////////////////////////////////////////////////////////////////

/* Update household first case and first infection time.
WARNING: this strucutre assumes a single first infector.
*/
void Household::update_first_case()
{
	double t_infec;
	double new_t0 = INFINITY; // NOTE: init to impossible values
	size_t new_case0 = -1;	  // NOTE: init to impossible values

	for (size_t i : infectors_list)
	{
		t_infec = (*indiv_ptr)[i].get_infec_time();
		new_t0 = std::min(t_infec, new_t0);
		if (t_infec == new_t0)
		{
			new_case0 = i;
		}
	}
	first_case = new_case0;
	start_followup = new_t0;
}

/* */ // FIXME : is then study_period a standard name? In general, how to load household data?
// FIXME : this should be a function of the time-models maybe
void Household::assign_end_followup()
{
	// find maximum study period within the household
	end_followup = 0;
	for (size_t i : infectors_list)
	{
		end_followup = std::max(end_followup, get_indiv(i).get_data("study_period"));
	}
}

/* Link an individual to the household */
void Household::link_individual(size_t i, int status)
{
	house_members.push_back(i);

	// status list
	if (status == HEALTHY)
	{
		susceptibles_list.push_back(i);
	}
	else if (status == SYMPTOMATIC)
	{
		infectors_list.push_back(i);
	}
	else
	{
		unknown_status_list.push_back(i);
	}
}

/* Change_status of individual i to new status (without controls). 
	WARNING : this introduce augmented variables -> it requires RJMCMC in case of inference.*/
void Household::change_status(size_t i, int new_status)
{
	int old_status = (*indiv_ptr)[i].get_status();

	if (old_status != new_status)
	{
		// erase individual i from previous list
		if (old_status == HEALTHY)
		{
			susceptibles_list.erase(find(susceptibles_list.begin(), susceptibles_list.end(), i));
		}
		else if (old_status == SYMPTOMATIC)
		{
			infectors_list.erase(find(infectors_list.begin(), infectors_list.end(), i));
		}
		else
		{
			unknown_status_list.erase(find(unknown_status_list.begin(), unknown_status_list.end(), i));
		}

		// add individual i to new list
		if (new_status == HEALTHY)
		{
			susceptibles_list.push_back(i);
		}
		else if (new_status == SYMPTOMATIC)
		{
			infectors_list.push_back(i);
		}
		else
		{
			unknown_status_list.push_back(i);
		}
	}
}