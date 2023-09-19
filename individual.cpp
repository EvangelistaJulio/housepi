// individual.cpp //

#include "individual.hpp"

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

void Base_Individual::add_augm(std::string augmID, Augmented new_augm)
{
	augm_var.emplace(augmID, new_augm);
}

void Base_Individual::add_augm(std::string augmID, Augmented new_augm, std::mt19937_64 &rng)
{
	augm_var.emplace(augmID, new_augm);
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
// //////////////// Household : (FIXME: temporary) /////////////////// //
////////////////////////////////////////////////////////////////////////

/* add individual */
void Household::add_individual(size_t i, int status)
{
	if (status == HEALTHY)
	{
		// Update the list of healthy residents within the household
		healthy_list.push_back(i);
	}
	else
	{
		// Update the list of symptomatic residents within the household
		symptomatic_list.push_back(i);
	}
}

/* change_status (UPDATE: too baroque) */
void Household::change_status(size_t i, int status)
{
	std::vector<size_t> *first;
	std::vector<size_t> *second;

	if (status == HEALTHY)
	{
		first = &healthy_list;
		second = &symptomatic_list;
	}
	else
	{
		first = &symptomatic_list;
		second = &healthy_list;
	}

	{
		// Find the index in the first vector
		auto it = std::find(first->begin(), first->end(), i);

		// if the index is in its old list erase it
		if (it != first->end())
		{
			// Erase the value from the first vector
			first->erase(it);
		}
	}

	{
		// Find the index in the first vector
		auto it = std::find(second->begin(), second->end(), i);

		// if the index is not the new list add it it
		if (it == second->end())
		{
			// Add the erased value to the second vector
			second->push_back(i);
		}
	}
}