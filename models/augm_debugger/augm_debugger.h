// augm_debugger.h //

#ifndef MODEL__H
#define MODEL__H

#include "mcmc.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// /////////////////////////////// Infectivity_Profile_Function //////////////////////////////// //
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Choice of the Infectivity Profile Function.
The infectivity profile is considered equal to the Time from Onset of Symptoms to Transmission `tost`
>	tost = <time infection infectee> - <time onset infector>
Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
where the pre-symptomaticvalues are rescaled as a function of the incubation time `t_incub`.
The arguments to provide are :
	0) "tost_alpha" : shape;
	1) "tost_mu" : location;
	2) "tost_sigma" : scale;
	3) "tost_tau" : mean incubation time;

Side note: symtoms can occur after infection (pre-symptomatic transmission), in other words `tost` can be negative.

----
Ref:
----
Ferretti, L. et al. The timing of COVID-19 transmission. 2020.09.04.20188516
Preprint at https://doi.org/10.1101/2020.09.04.20188516 (2020).*/

/* Define the Infectivity_Profile_Function class */
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// ////////////////////////////////////////// Model : class /////////////////////////////////////////// //
//////////////////////////////////////////////////////////////////////////////////////////////////////////

class Model : public MCMC
{
protected:
	std::vector<double> inst_profile_LUT = {};
	std::vector<double> cumul_profile_LUT = {};

public:
	Model() : infec_profile_func(){};

	//=====Model Specifc Methods =====//
	Infectivity_Profile_Function infec_profile_func;
	double compute_infec_profile_inst(int i, double t_transm);
	double compute_infec_profile_cumul(int i, double t_transm);
	double hazards_of_infec(std::vector<double> *profiles, int idx);

	void augm_var_store_dependecies(int i) override;
	double get_infec_time(int i) const override;
	double compute_single_house_loglike(int h) override;
	void initialize_model() override;

	void simulate_single_house(int h) override;
	void format_augmented_data() override;

};

#endif // MODEL__H