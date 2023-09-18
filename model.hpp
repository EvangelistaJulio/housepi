// model.hpp //

#ifndef MODEL__H
#define MODEL__H

#include <string>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <random>
#include <iostream>

#include "household.hpp"

///////////////////////////////////////////////////////////////////////
// //////////////// Onset to Transmission : class ////////////////// //
///////////////////////////////////////////////////////////////////////

class Onset_to_Transm : public Base_Distrib
{
public:
    // Constructors
    Onset_to_Transm(kwargs new_args) : Base_Distrib(new_args) {}
    Onset_to_Transm(kwargs_ext new_args_ext) : Base_Distrib(new_args_ext) {}
    ~Onset_to_Transm(){};

    double pdf(timetype tost, timetype t_incub);
    double cdf(timetype tost, timetype t_incub);
};


///////////////////////////////////////////////////////////////////////
// ///////////////////////// Model : class ///////////////////////// //
///////////////////////////////////////////////////////////////////////

class Model : public Parameter_Dictionary
{
private:
	Onset_to_Transm * infectivity_profile;

public:
	Model() {}
	~Model() {}

	// Augmetend Data 
	void assign_augm_data(Household *H, std::mt19937_64 &rng);
	void _assign_augm_to_ind_(Household *H, size_t idx, std::mt19937_64 &rng); // UPDATE: temporary
	timetype get_infec_time(Household *H, size_t idx) const;

	// Infectivity Profile
	void assign_infec_profile();
	double infectivity_inst(Household *H, size_t idx, timetype t_transm);
	double infectivity_cumul(Household *H, size_t idx, timetype t_transm);
	void init_infectivity_profiles(Household *H);
	void update_infectivity_profiles(Household *H, size_t idx);
	void update_first_case(Household *H);

	double risk_beta(Household *H, const std::vector<std::vector<double>> &rate_vector, size_t idx);
	double log_surv_proba(Household *H, size_t idx);
	double log_infec_proba(Household *H, size_t idx);
	double compute_loglike(Household *H);

	double community_risk_inst(timetype t1);
	double community_risk_cumul(timetype t0, timetype t1);
	double contact_modulation(Household *H, size_t infectee, size_t idx);
	double susceptibility(Household *H, size_t idx);
	double transmissability(Household *H);

	void simulate_data(Household *H, timetype delta_t, std::mt19937_64 &rng);
};

#endif // MODEL__H