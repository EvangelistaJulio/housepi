/*
 * profiles.cpp
 *
 * Copyright (C) 2024 Francesco Camaglia, Institut Pasteur
 */

#include "profiles.h"
#include "distribution.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////// Base_Infectivity_Profile : class ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

double Base_Infectivity_Profile::pdf(double serial, const std::vector<double> &x) const
{
	UNUSED(serial);
	UNUSED(x);
	tools::error("Call to undefined method 'pdf' of infectivity profile '" + get_name() + "'.");
	return std::numeric_limits<double>::quiet_NaN();
}

double Base_Infectivity_Profile::cdf(double serial, const std::vector<double> &x) const
{
	UNUSED(serial);
	UNUSED(x);
	tools::error("Call to undefined method 'pdf' of infectivity profile '" + get_name() + "'.");
	return std::numeric_limits<double>::quiet_NaN();
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Infec_Prof_Ferretti20 : class //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
The infectivity profile is considered equal to the distribution
of the Time from Onset of Symptoms to Transmission 't_ost':
>	t_ost = <time infection infectee> - <time onset infector>
Its distribution is modeled as a skew-logistic distribution (generalized logistic type I),
where the pre-symptomatic values are rescaled as a function of the incubation time 'incub'.
The arguments to provide are :
  0) "alpha" : shape;
  1) "mu" : location;
  2) "sigma" : scale;
  3) "tau" : mean incubation time;
Note: symtoms can occur after infection (pre-symptomatic transmission),
i. e. 't_ost' can be negative.
----
Ref:
----
Ferretti L, Ledda A, Wymant C, Zhao L, Ledda V, Abeler-Dorner L, et al.
The timing of COVID-19 transmission.
medRxiv; 2020. p. 2020.09.04.20188516. doi:10.1101/2020.09.04.20188516
 */
class Infec_Prof_Ferretti20 : public Base_Infectivity_Profile
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

	//-----  Probability distribution function -----//
	double pdf(double serial, const std::vector<double> &x) const override
	{
		double incub = x[0];		   // incubation period
		double t_OST = serial - incub; // t. from onset to transmission

		double a = get_arg(0); // shape argument
		double m = get_arg(1); // location argument
		double s = get_arg(2); // scale argument
		double T = get_arg(3); // mean incubation period argument

		double norm_pdf = (s / a) * (1. + (incub / T - 1.) * pow(1 + exp(m / s), -a) - (incub / T) * (pow(1. + exp((m + T) / s), -a)));

		if (t_OST < 0)
		{
			t_OST *= (T / incub);
		}
		return exp((m - t_OST) / s) * pow(1. + exp((m - t_OST) / s), -1. - a) / norm_pdf;
	}

	//-----  Cumulative distribution function -----//
	double cdf(double serial, const std::vector<double> &x) const override
	{
		double incub = x[0];		   // incubation period
		double t_OST = serial - incub; // t. from onset to transmission

		double a = get_arg(0); // shape argument
		double m = get_arg(1); // location argument
		double s = get_arg(2); // scale argument
		double T = get_arg(3); // mean incubation period argument

		double norm_cdf = (1. + (incub / T - 1.) * pow(1 + exp(m / s), -a) - (incub / T) * (pow(1. + exp((m + T) / s), -a)));

		if (t_OST < 0)
		{
			// Integrate from -incub to t_ost
			return (incub / T) * (pow(1. + exp((m - t_OST * T / incub) / s), -a) - pow(1. + exp((m + T) / s), -a)) / norm_cdf;
		}
		else
		{
			// Integrate from -incub to 0
			double tmp = (incub / T) * (pow(1. + exp(m / s), -a) - pow(1. + exp((m + T) / s), -a));
			// Integrate from 0 to t_ost
			return (tmp + pow(1. + exp((m - t_OST) / s), -a) - pow(1. + exp(m / s), -a)) / norm_cdf;
		}
	}
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Infec_Prof_Ashcroft20 : class //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
 ----
 Ref:
 ----
 Ashcroft P, Huisman JS, Lehtinen S, Bouman JA, Althaus CL, Regoes RR, et al.
 COVID-19 infectivity profile correction. Swiss medical weekly. 2020. p. w20336.
*/
class Infec_Prof_Ashcroft20 : public Base_Infectivity_Profile
{
private:
	std::shared_ptr<Base_Distribution> gamma_distr_ptr;

public:
	// Constructor with immutable values
	Infec_Prof_Ashcroft20(const dictargs &new_args = {})
	{
		init("ashcroft20", {{"shape", 97.32}, {"scale", 0.268}, {"shift", 25.6}},
			 {"shape", "scale", "shift"});
		assign_args(new_args);
		dictargs aux_args_val = {{"shape", get_arg(0)}, {"scale", get_arg(1)}};
		gamma_distr_ptr = distribution("gamma", aux_args_val, {{"min", get_arg(2)}});
	}
	// Deconstructor
	~Infec_Prof_Ashcroft20() {};

	double pdf(double serial, const std::vector<double> &x) const override
	{
		double incub = x[0];		   // incubation period
		double t_OST = serial - incub; // t. from onset to transmission

		double shift = get_arg(2); // shift argument
		return gamma_distr_ptr->pdf(shift + t_OST);
	}

	double cdf(double serial, const std::vector<double> &x) const override
	{
		double incub = x[0];		   // incubation period
		double t_OST = serial - incub; // t. from onset to transmission

		double shift = get_arg(2); // shift argument
		return gamma_distr_ptr->cdf(shift + t_OST);
	}
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

//* Return a pointer to infectivity profile */
std::unique_ptr<Base_Infectivity_Profile> infectivity_profile_switchboard(const std::string &name)
{
	std::unique_ptr<Base_Infectivity_Profile> output = std::make_unique<Base_Infectivity_Profile>();
	const std::vector<std::string> options = {"ferretti20", "ashcroft20"};

	if (name == "ferretti20")
	{
		output = std::make_unique<Infec_Prof_Ferretti20>();
	}
	else if (name == "ashcroft20")
	{
		output = std::make_unique<Infec_Prof_Ashcroft20>();
	}
	else if (!tools::is_na(name))
	{
		std::string err_message = "Unknown infectivity profile '" + name + "'.";
		err_message += tools::find_closest(name, options);
		throw std::invalid_argument(err_message);
	}
	return output;
}