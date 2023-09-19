// individual.hpp //

#ifndef INDIVIDUAL__H
#define INDIVIDUAL__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "varia.hpp"

// model
#define HEALTHY 0
#define SYMPTOMATIC 1

/////////////////////////////////////////////////////////
///////// Base_Individual : Class declaration ///////////////
////////////////////////////////////////////////////////

class Base_Individual
{
private:
	//=====Attributes=====//
	int status = 0;
	int index = 0;

	std::unordered_map<std::string, double> data = {};
	std::unordered_map<std::string, Augmented> augm_var = {};

public:
	//=====Constructors=====//
	Base_Individual(
		int status,
		int index,
		std::unordered_map<std::string, double> data);
	Base_Individual() = default;
	~Base_Individual(){};

	//=====Methods=====//

	//=====Getter=====//
	bool get_status() const { return status; }
	int get_index() const { return index; }
	double get(std::string dataID) const { return data.at(dataID); }
	double get_augm(std::string augmID) const { return augm_var.at(augmID).get_value(); }
	bool augm_infer(std::string augmID) { return augm_var.at(augmID).get_infer(); }

	//=====Setter=====//
	void add_augm(std::string augmID, Augmented new_augm, bool init = false);

	void set_status(int status) { this->status = status; }
	void set_augm(std::string augmID, double new_val) { augm_var.at(augmID).set_value(new_val); }
	void init_augm(std::string augmID, std::mt19937_64 &rng) { set_augm(augmID, augm_model_rand(augmID, rng)); }

	//=====model distribution=====//
	double augm_model_rand(std::string augmID, std::mt19937_64 &rng) { return augm_var.at(augmID).model_rand(rng); }
	double augm_model_log_pdf(std::string augmID);
	double augm_model_log_ratio(std::string augmID, double old_value, double new_val);

	//=====proposer=====//
	double augm_propose(std::string augmID, std::mt19937_64 &rng);
	double augm_propose_log_ratio(std::string augmID, double old_value, double new_val);
	bool augm_valid(std::string augmID, double new_val) { return augm_var.at(augmID).valid(new_val); }
};

/////////////////////////////////////////////////////////
///////// Household : Class declaration ///////////////
////////////////////////////////////////////////////////

class Household
{
public:
	// Attributes
	int house_index = -1;
	size_t size = 0;
	size_t first_case = -1;
	timetype infec_time = INFINITY;
	std::vector<size_t> healthy_list = {};
	std::vector<size_t> symptomatic_list = {};
	double loglike = 0;

	// Constructors
	Household() {}
	~Household(){};

	// Getter
	int get_index() const { return house_index; }

	// Setter
	void set_index(int index) { house_index = index; }
	void add_individual(size_t i, int status);
	void change_status(size_t i, int status);
};

#endif // INDIVIDUAL__H