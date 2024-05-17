// single_introduction_augmented.h //
// Copyright (C) 2024 Francesco Camaglia, Institut Pasteur //

/* 
//////////////////////////////////////////////////////////////////////////////////////////////////////////

Choice of the Infectivity Profile Function.
-----
The infectivity profile is considered equal to the Time from Onset of Symptoms to Transmission `tost`
>	tost = <time infection infectee> - <time onset infector>
Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
where the pre-symptomatic values are rescaled as a function of the incubation time `t_incub`.
The arguments to provide are :
	0) "tost_alpha" : shape;
	1) "tost_mu" : location;
	2) "tost_sigma" : scale;
	3) "tost_tau" : mean incubation time;
Side note: symtoms can occur after infection (pre-symptomatic transmission), i. e. `tost` can be negative.
----
Ref:
----
Ferretti L, Ledda A, Wymant C, Zhao L, Ledda V, Abeler-Dörner L, et al. 
The timing of COVID-19 transmission. 
medRxiv; 2020. p. 2020.09.04.20188516. doi:10.1101/2020.09.04.20188516
//////////////////////////////////////////////////////////////////////////////////////////////////////////
*/

#ifndef MODEL__H
#define MODEL__H

#include "mcmc.h"

// test with the epidemic curve
//#define TEST_EPI_CURVE

// alias for a square matrix
typedef std::vector<std::vector<double>> squarematrix; 

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Infectivity_Profile_Function ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Define the Infectivity_Profile_Function class. */
class Infectivity_Profile_Function : public Base_Function
{
public:
	// Constructors
	Infectivity_Profile_Function() {};
	Infectivity_Profile_Function(sdargs new_args) : Base_Function(new_args, {"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}) {}
	Infectivity_Profile_Function(sdargs_ext new_args_ext) : Base_Function(new_args_ext, {"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}) {}
	~Infectivity_Profile_Function(){};

	double pdf(double tost, double t_incub);
	double cdf(double tost, double t_incub);
};

//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------//

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
public:
	//=====Attributes for Look-Up-Tables=====//
	std::vector<squarematrix> inst_profile_LUT = {};
	std::vector<squarematrix> cumul_profile_LUT = {};

	//=====Constructor=====//
	Model() : infec_profile_func(){};

#ifdef TEST_EPI_CURVE
	void build_extra(std::string extra_pathtofile ) override;
#endif

	//=====Infectivity Profiles=====//
	Infectivity_Profile_Function infec_profile_func;
	double compute_infec_profile_inst(int h, int i, double t_transm);
	double compute_infec_profile_cumul(int h, int i, double t_transm);
	double hazards_of_infec(int h, int i, std::vector<squarematrix> &profiles);

	//=====Override General Methods=====//
	void initialize_model() override;
	void augm_var_store_dependecies(int h, int i) override;
	double get_infec_time(int h, int i) const override;

	//=====Override Methods for Inference=====//
	double compute_single_house_loglike(int h) override;

	//=====Rate Choices=====//
	void assign_rate_per_group_to_individual(Individual &this_indiv) override;
	void assign_rate_per_contact_to_household(int h, int from, int to) override;

	//=====Override Methods for Simulations=====//
	void simulate_single_house(int h) override;
	void format_augmented_data() override;
};

#endif // MODEL__H