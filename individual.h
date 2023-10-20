// individual.h //

#ifndef INDIVIDUAL__H
#define INDIVIDUAL__H

#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "varia.h"

// model
#define HEALTHY 0
#define SYMPTOMATIC 1

/////////////////////////////////////////////////////////
///////// Individual : Class declaration ///////////////
////////////////////////////////////////////////////////

class Individual
{
private:
	//=====Attributes=====//
	int status = 0;
	int index = 0;

	std::unordered_map<std::string, double> data = {};
	std::unordered_map<std::string, Augmented> augm_var = {};

	const double *_infecti = nullptr;
	const double *_suscept = nullptr;

public:
	size_t home = -1;

	//=====Constructors=====//
	Individual(int status, int index, std::unordered_map<std::string, double> data);
	Individual() = default;
	~Individual(){};

	//=====Methods for the Model=====//
	timetype get_infec_time() const;

	//=====Getter=====//
	bool get_status() const { return status; }
	int get_index() const { return index; }
	double get_data(std::string dataID) const
	{
		// UPDATE : this could be a bottleneck if data are called every iteration
		if (data.find(dataID) == data.end())
		{
			std::cerr << "\nERROR: data column `" << dataID << "` unknown." << std::endl;
			std::cerr << "Exit...\n"
					  << std::endl;
			std::exit(EXIT_FAILURE);
		}
		return data.at(dataID);
	}
	double get_augm(std::string augmID) const { return augm_var.at(augmID).get_value(); }
	bool augm_infer(std::string augmID) { return augm_var.at(augmID).get_infer(); }

	double get_suscept() const { return *_suscept; }
	double get_infecti() const { return *_infecti; }

	//=====Setter=====//
	void add_augm(std::string augmID, Augmented new_augm);
	void add_augm(std::string augmID, Augmented new_augm, std::mt19937_64 &rng);

	void set_status(int d) { status = d; }
	void set_augm(std::string augmID, double new_val) { augm_var.at(augmID).set_value(new_val); }
	void init_augm(std::string augmID, std::mt19937_64 &rng);

	void set_suscept(const double *x) { _suscept = x; }
	void set_infecti(const double *x) { _infecti = x; }

	//=====Model distribution=====//
	double augm_model_rand(std::string augmID, std::mt19937_64 &rng)
	{
		return augm_var.at(augmID).model_rand(rng);
	}
	double augm_model_log_pdf(std::string augmID);
	double augm_model_log_ratio(std::string augmID, double old_value, double new_val);

	//=====Proposer=====//
	double augm_propose(std::string augmID, std::mt19937_64 &rng);
	double augm_propose_log_ratio(std::string augmID, double old_value, double new_val);
	bool augm_valid(std::string augmID, double new_val) { return augm_var.at(augmID).valid(new_val); }
};

/////////////////////////////////////////////////////////
////////// Household : Class declaration ///////////////
////////////////////////////////////////////////////////

class Household
{
private:
	Individual **indiv_ptr = nullptr;
	double loglike = -INFINITY;
	double buffer = -INFINITY;

public:
	//=====Attributes=====//
	int house_index = -1;
	size_t size = 0;
	size_t first_case = -1;
	timetype infec_time = INFINITY;
	std::vector<size_t> susceptibles_list = {};
	std::vector<size_t> infectors_list = {};

	//=====Constructors=====//
	Household() {}
	~Household(){};

	//======Getter======//
	int get_index() const { return house_index; }

	//======Setter======//
	void set_index(int index) { house_index = index; }
	void set_indiv_ptr(Individual **new_ptr) { indiv_ptr = new_ptr; }
	void add_individual(size_t i, int status);
	void change_status(size_t i, int new_status);
	void update_first_case();

	//======Loglike======//
	double get_loglike() const { return loglike; }
	void set_loglike(double new_loglike)
	{
		buffer = loglike;
		loglike = new_loglike;
	}
	void refuse() { loglike = buffer; }
};

#endif // INDIVIDUAL__H