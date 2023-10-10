// individual.cpp //

#include "individual.h"

using namespace std;

////////////////////////////////////////////////////////////////////////
// //////////////// Base_Individual : Constructor /////////////////// //
////////////////////////////////////////////////////////////////////////

Base_Individual::Base_Individual(
	int status,
	int index,
	std::unordered_map<std::string, double> data)
{
	this->status = status;
	this->index = index;
	this->data = data;
}

/* Add an augmented variable to the individual */
void Base_Individual::add_augm(std::string augmID, Augmented new_augm)
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
	}
}

/* Init augmented variable within a random initial value */
void Base_Individual::init_augm(std::string augmID, std::mt19937_64 &rng)
{
	double init_value = augm_model_rand(augmID, rng);
	while (!augm_valid(augmID, init_value))
	{
		init_value = augm_model_rand(augmID, rng);
	}
	set_augm(augmID, init_value);
}

/* Add an augmented variable to the individual, plus initialize it */
void Base_Individual::add_augm(std::string augmID, Augmented new_augm, std::mt19937_64 &rng)
{
	add_augm(augmID, new_augm);
	init_augm(augmID, rng);
}

double Base_Individual::augm_model_log_pdf(std::string augmID)
{
	double new_val = augm_var.at(augmID).get_value();
	return augm_var.at(augmID).model_log_pdf(new_val);
}

double Base_Individual::augm_model_log_ratio(std::string augmID, double old_val, double new_val)
{
	return augm_var.at(augmID).model_log_ratio(old_val, new_val);
}

double Base_Individual::augm_propose(std::string augmID, std::mt19937_64 &rng)
{
	return augm_var.at(augmID).propose(rng);
}

double Base_Individual::augm_propose_log_ratio(std::string augmID, double old_val, double new_val)
{
	return augm_var.at(augmID).propose_log_ratio(old_val, new_val);
}

////////////////////////////////////////////////////////////////////////
// /////////////////////////// Household //////////////////////////// //
////////////////////////////////////////////////////////////////////////

// FIXME: this architecture fixes the possible status of the individual

/* add individual */
void Household::add_individual(size_t i, int status)
{
	(status == HEALTHY) ? healthy_list.push_back(i) : symptomatic_list.push_back(i);
}

/* change_status of individual i to new status (without controls) */
void Household::change_status(size_t i, int new_status)
{
	std::vector<size_t> *new_list = (new_status == HEALTHY) ? &healthy_list : &symptomatic_list;
	std::vector<size_t> *old_list = (new_status == HEALTHY) ? &symptomatic_list : &healthy_list;

	// erase the index from the first and add to the second
	vector<size_t>::iterator element = std::find(old_list->begin(), old_list->end(), i);

	old_list->erase(element);
	new_list->push_back(i);
}