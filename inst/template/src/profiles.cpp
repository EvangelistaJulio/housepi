/*
 * profiles.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "profiles.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Infec_Prof_Ferretti20 : probability distribution  ////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Probability distribution function for the Infectivity Profile. */
double Infec_Prof_Ferretti20::pdf(double tost, double t_incub) const
{		
	double a = get_arg(0); // shape argument
	double m = get_arg(1); // location argument
	double s = get_arg(2); // scale argument
	double T = get_arg(3); // mean incubation period argument

	double norm_pdf = (s / a) * (1. + (t_incub / T - 1.) * pow(1 + exp(m / s), -a) - (t_incub / T) * (pow(1. + exp((m + T) / s), -a)));

	if (tost < 0)
	{
		tost *= (T / t_incub);
	}
	return exp((m - tost) / s) * pow(1. + exp((m - tost) / s), -1. - a) / norm_pdf;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////// Infec_Prof_Ferretti20 : cumulative distribution  ////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Cumulative distribution function forthe Infectivity Profile. */
double Infec_Prof_Ferretti20::cdf(double tost, double t_incub) const
{
	double a = get_arg(0); // shape argument
	double m = get_arg(1); // location argument
	double s = get_arg(2); // scale argument
	double T = get_arg(3); // mean incubation period argument

	double norm_cdf = (1. + (t_incub / T - 1.) * pow(1 + exp(m / s), -a) - (t_incub / T) * (pow(1. + exp((m + T) / s), -a)));

	if (tost < 0)
	{
		// Integrate from -t_incub to tost
		return (t_incub / T) * (pow(1. + exp((m - tost * T / t_incub) / s), -a) - pow(1. + exp((m + T) / s), -a)) / norm_cdf;
	}
	else
	{
		// Integrate from -t_incub to 0
		double tmp = (t_incub / T) * (pow(1. + exp(m / s), -a) - pow(1. + exp((m + T) / s), -a));
		// Integrate from 0 to tost
		return (tmp + pow(1. + exp((m - tost) / s), -a) - pow(1. + exp(m / s), -a)) / norm_cdf;
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////// Infec_Prof_Ashcroft20 : probability distribution  ////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

/* Probability distribution function for the Infectivity Profile. */
double Infec_Prof_Ashcroft20::pdf(double t_infection) const
{
	double shift = get_arg(2); // shift argument

	return gamma_distr_ptr->pdf(shift + t_infection);
}

/* Cumulative distribution function for the Infectivity Profile. */
double Infec_Prof_Ashcroft20::cdf(double t_infection) const
{
	double shift = get_arg(2); // shift argument

	return gamma_distr_ptr->cdf(shift + t_infection);
}
