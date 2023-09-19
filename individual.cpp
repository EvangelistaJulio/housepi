// individual.cpp //

#include "individual.hpp"

using namespace std;

/////////////////////////////////////////////////////
//  Base_Individual : Class
/////////////////////////////////////////////////////

// ----------------------------
//  Constructors
// ----------------------------

Base_Individual::Base_Individual(
	int status,
	int index,
	std::unordered_map<std::string, double> data)
{
	this->status = status;
	this->index = index;
	this->data = data;
}

void Base_Individual::add_augm(std::string augmID,
							   double new_val,
							   std::string model_name,
							   kwargs_ext model_args_ext,
							   std::string strategy_name,
							   kwargs strategy_args,
							   double value_0,
							   kwargs value_range)
{
	Augmented this_augm(new_val, model_name, model_args_ext, strategy_name, strategy_args, value_0, value_range);
	augm_var.emplace(augmID, this_augm);
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
