/*
 * profiles.h
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#ifndef PROFILES__H
#define PROFILES__H

#include "varia.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Infec_Prof_Ferretti20 : class ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
The infectivity profile is considered equal to the distribution
of the Time from Onset of Symptoms to Transmission 'tost':
>	tost = <time infection infectee> - <time onset infector>
Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
where the pre-symptomatic values are rescaled as a function of the incubation time 't_incub'.
The arguments to provide are :
  0) "alpha" : shape;
  1) "mu" : location;
  2) "sigma" : scale;
  3) "tau" : mean incubation time;
Note: symtoms can occur after infection (pre-symptomatic transmission),
i. e. 'tost' can be negative.
----
Ref:
----
Ferretti L, Ledda A, Wymant C, Zhao L, Ledda V, Abeler-Dorner L, et al.
The timing of COVID-19 transmission.
medRxiv; 2020. p. 2020.09.04.20188516. doi:10.1101/2020.09.04.20188516
 */
class Infec_Prof_Ferretti20 : public Base_Function
{
public:
    // Constructor with immutable values
    Infec_Prof_Ferretti20(const dictargs &new_args = {})
    {
        init("ferretti20", {{"alpha", 5.85}, {"mu", -4.0}, {"sigma", 1.85}, {"tau", 5.8}},
             {"alpha", "mu", "sigma", "tau"});
        assign_args(new_args);
    }
    // Deconstructor
    ~Infec_Prof_Ferretti20() = default;

    double pdf(double t_onset_to_transmission, double incubation_period) const;
    double cdf(double t_onset_to_transmission, double incubation_period) const;
};

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Infec_Prof_Ashcroft20 : class ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
 ----
 Ref:
 ----
 Ashcroft P, Huisman JS, Lehtinen S, Bouman JA, Althaus CL, Regoes RR, et al.
 COVID-19 infectivity profile correction. Swiss medical weekly. 2020. p. w20336.
*/
class Infec_Prof_Ashcroft20 : public Base_Function
{
private:
    Base_Distribution *gamma_distr_ptr;

public:
    // Constructor with immutable values
    Infec_Prof_Ashcroft20(const dictargs &new_args = {})
    {
        init("ashcroft20", {{"shape", 97.32}, {"scale", 0.268}, {"shift", 25.6}},
             {"shape", "scale", "shift"});
        assign_args(new_args);
        dictargs aux_args_val = {{"shape", get_arg(0)}, {"scale", get_arg(1)}};
        gamma_distr_ptr = distribution_ptr("gamma", aux_args_val, {{"min", get_arg(2)}});
    }
    // Deconstructor
    ~Infec_Prof_Ashcroft20() {};

    double pdf(double t_infection) const;
    double cdf(double t_infection) const;
};

#endif // PROFILES__H