/*
 * infectivity_profiles.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef PROFILES__H
#define PROFILES__H

#include "utils.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Infec_Prof_Skew_Logistic : class ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*

The infectivity profile is considered equal to the distribution
of the Time from Onset of Symptoms to Transmission `tost`:
>	tost = <time infection infectee> - <time onset infector>
Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
where the pre-symptomatic values are rescaled as a function of the incubation time `t_incub`.
The arguments to provide are :
	0) "tost_alpha" : shape;
	1) "tost_mu" : location;
	2) "tost_sigma" : scale;
	3) "tost_tau" : mean incubation time;
Note: symtoms can occur after infection (pre-symptomatic transmission),
i. e. `tost` can be negative.
----
Ref:
----
Ferretti L, Ledda A, Wymant C, Zhao L, Ledda V, Abeler-Dörner L, et al.
The timing of COVID-19 transmission.
medRxiv; 2020. p. 2020.09.04.20188516. doi:10.1101/2020.09.04.20188516
///////////////////////////////////////////////////////////////////////////////////////////////////
*/

/* Define the Inf_Prof_Skew_Logistic class. */
class Infec_Prof_Skew_Logistic : public Base_Function
{
public:
	// Constructors
	Infec_Prof_Skew_Logistic() {};
	Infec_Prof_Skew_Logistic(sdargs new_args) : Base_Function(new_args, {"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}) {}
	Infec_Prof_Skew_Logistic(sdargs_ptr new_args_ext) : Base_Function(new_args_ext, {"tost_alpha", "tost_mu", "tost_sigma", "tost_tau"}) {}
	~Infec_Prof_Skew_Logistic() {};

	double pdf(double tost, double t_incub);
	double cdf(double tost, double t_incub);
};

#endif // PROFILES__H