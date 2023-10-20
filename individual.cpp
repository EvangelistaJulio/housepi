// individual.cpp //

#include "individual.h"

using namespace std;

#define MAX_NUMBER_ATTEMPTS 10

////////////////////////////////////////////////////////////////////////
// //////////////// Individual : Constructor /////////////////// //
////////////////////////////////////////////////////////////////////////

Individual::Individual(
	int new_status,
	int new_index,
	std::unordered_map<std::string, double> new_data)
{
	status = new_status;
	index = new_index;
	data = new_data;
}

/* Init augmented variable within a random initial value */
void Individual::init_augm(std::string augmID, std::mt19937_64 &rng)
{
	int num_attempts = 0;
	double init_value = augm_model_rand(augmID, rng);

	while ((!augm_valid(augmID, init_value)))
	{
		if (num_attempts >= MAX_NUMBER_ATTEMPTS)
		{
			std::cerr << "ERROR: Number of attempts exceed in random initialization of " << augmID << std::endl;
			std::cerr << "Exit...\n"
					  << std::endl;
			exit(EXIT_FAILURE);
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
		std::cerr << "ERROR: Augmented variable `" << augmID << "` has multiple definitions." << std::endl;
		std::cerr << "Exit...\n"
				  << std::endl;
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

/* Update household first case and first infection time.
WARNING: this strucutre assumes a single first infector.
*/
void Household::update_first_case()
{
	timetype t_infec;
	timetype new_t0 = INFINITY;
	size_t new_case0 = INFINITY; // NOTE: init to impossible values

	for (size_t i : infectors_list)
	{
		t_infec = (*indiv_ptr)[i].get_infec_time();
		new_t0 = std::min(t_infec, new_t0);
		if (t_infec == new_t0)
		{
			new_case0 = i;
		}
	}
	infec_time = new_t0;
	first_case = new_case0;
}

/* Add individual to the household */
void Household::add_individual(size_t i, int status)
{
	(status == HEALTHY) ? susceptibles_list.push_back(i) : infectors_list.push_back(i);
}

/* Change_status of individual i to new status (without controls). 
WARNING : this introduce augmented variables -> it requires RJMCMC in case of inference.*/
void Household::change_status(size_t i, int new_status)
{
	std::vector<size_t> *new_list = (new_status == HEALTHY) ? &susceptibles_list : &infectors_list;
	std::vector<size_t> *old_list = (new_status == HEALTHY) ? &infectors_list : &susceptibles_list;

	// erase the index from the first and add to the second
	vector<size_t>::iterator element = std::find(old_list->begin(), old_list->end(), i);

	old_list->erase(element);
	new_list->push_back(i);
}