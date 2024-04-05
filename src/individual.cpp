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
	/* // FIXME : control should be at the point of pre-processing
	if (data.find(dataID) == data.end())
	{
		std::string err_text = "Data column `" + dataID + "` unknown.";
		raise_error_message(err_text);
		std::exit(EXIT_FAILURE);
	}
	*/
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
		raise_error_message(err_text);
		std::exit(EXIT_FAILURE);
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
			std::string err_text = "Number of attempts exceeded in random initialization of `" + augmID + "`.";
			raise_error_message(err_text);
			std::exit(EXIT_FAILURE);
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
		raise_error_message(err_text);
		std::exit(EXIT_FAILURE);
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

/* Assign a data feature to the household */
void Household::assign_data(std::string dataID, double value)
{
	if (house_data.find(dataID) == house_data.end())
	{
		// if rateID was undeclared, add it
		house_data.emplace(dataID, value);
	}
	else
	{
		// rateID was previously declared
		std::string err_text = "Found household with multiple definitions for data `" + dataID + "`.";
		raise_error_message(err_text);
		std::exit(EXIT_FAILURE);
	}
}

/* get_data : bottleneck if data are called every iteration. Deprecated outisde of pre-processing.*/
double Household::get_data(std::string dataID) const
{
	/* // FIXME : control should be done at the moment of pre-processing
	if (house_data.find(dataID) == house_data.end())
	{
		std::string err_text = "Household data  `" + dataID + "` unknown.";
		raise_error_message(err_text);
		std::exit(EXIT_FAILURE);
	}
	*/
	return house_data.at(dataID);
}

/* Link an individual to the household */
void Household::link_individual(int i, int status)
{
	// status list
	if (status == SYMPTOMATIC)
	{
		infectious_list.push_back(i);
	}
	else if (status == HEALTHY)
	{
		susceptibles_list.push_back(i);
	}
	else
	{
		unknown_status_list.push_back(i);
	}
	house_size++;
}

/* Return the indexes of the infectious, then the susceptibles, unknown statuses */
std::vector<int> Household::get_house_members()
{
	std::vector<int> members;
	members.reserve(house_size); // Reserve space

	// Copy infectious
	members.insert(members.end(), infectious_list.begin(), infectious_list.end());

	// Copy susceptibles
	members.insert(members.end(), susceptibles_list.begin(), susceptibles_list.end());

	// Copy unknown
	members.insert(members.end(), unknown_status_list.begin(), unknown_status_list.end());

	return members;
}

/* Change_status of individual i to new status (without controls).
	WARNING : this introduce augmented variables -> it requires RJMCMC in case of inference.*/
void Household::change_status(int i, int new_status)
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
			infectious_list.erase(find(infectious_list.begin(), infectious_list.end(), i));
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
			infectious_list.push_back(i);
		}
		else
		{
			unknown_status_list.push_back(i);
		}

		// change the status of individual to new
		(*indiv_ptr)[i].set_status(new_status);
	}
}

// update the values for the Bayesian diagnostics
void Household::update_dgn()
{
	dgn_mean_ll = (loglike + dgn_mean_ll * dgn_sample_size) / (dgn_sample_size + 1);
	dgn_mean_sqr_ll = (std::pow(loglike, 2) + dgn_mean_sqr_ll * dgn_sample_size) / (dgn_sample_size + 1);
	dgn_sample_size++;
}