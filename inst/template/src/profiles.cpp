/*
 * profiles.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "mcmc.h"

Base_InfecProfile::Base_InfecProfile() {}
Base_InfecProfile::Base_InfecProfile(std::function<double(parID)> get) : _external_param_getter_(std::move(get)) {}
Base_InfecProfile::~Base_InfecProfile() = default;

void Base_InfecProfile::init(const std::string &name, const dict_args &defaults)
{
	this->_name_ = name;
	for (const auto &pair : defaults)
	{
		this->_arg_names_.emplace_back(pair.first);
		this->_values_.emplace_back(pair.second);
	}
}

double Base_InfecProfile::pdf(const std::vector<double> &x)
{
	UNUSED(x);
	throw std::invalid_argument(
		"Call to undefined method 'pdf' of infectivity profile '" + _name_ + "'.");
}
double Base_InfecProfile::cdf(const std::vector<double> &x)
{
	UNUSED(x);
	throw std::invalid_argument(
		"Call to undefined method 'cdf' of infectivity profile '" + _name_ + "'.");
}


/* Return value indexed as 'i'. */
double Base_InfecProfile::get_arg(ushort i)
{
	return (this->*_func_get_arg_)(i);
}

/* Return internal value indexed as 'i'. */
double Base_InfecProfile::_internal_get_arg_(ushort i)
{
	return _values_[i];
}

/* Return external value indexed as 'i' with functional . */
double Base_InfecProfile::_external_get_arg_(ushort i)
{
	return _external_param_getter_(_indexes_[i]);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// InfecProfile_Ferretti20 : class //////////////////////////////////
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
class InfecProfile_Ferretti20 : public Base_InfecProfile
{
private:
	void defaults()
	{
		init("ferretti20", {{"alpha", 5.85}, {"mu", -4.0}, {"sigma", 1.85}, {"tau", 5.8}});
	}

public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Ferretti20() : Base_InfecProfile() { defaults(); };
	InfecProfile_Ferretti20(std::function<double(parID)> get) : Base_InfecProfile(get) { defaults(); }

	//*  Probability distribution function. *//
	double pdf(const std::vector<double> &x) override
	{
		double t_OST = x[0] - x[1]; // t. from onset to transmission
		double incub = x[1];		// incubation period

		double a = get_arg(0); // shape argument
		double m = get_arg(1); // location argument
		double s = get_arg(2); // scale argument
		double T = get_arg(3); // mean incubation period argument

		double norm_pdf = (s / a) * (1. + (incub / T - 1.) * pow(1 + exp(m / s), -a) - (incub / T) * pow(1. + exp((m + T) / s), -a));

		if (t_OST < 0)
		{
			t_OST *= (T / incub);
		}

		return exp((m - t_OST) / s) * pow(1. + exp((m - t_OST) / s), -a - 1.) / norm_pdf;
	}

	//*  Cumulative distribution function. *//
	double cdf(const std::vector<double> &x) override
	{
		double t_OST = x[0] - x[1]; // t. from onset to transmission
		double incub = x[1];		// incubation period

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
////////////////////////////////// InfecProfile_Ashcroft20 : class //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*
 ----
 Ref:
 ----
 Ashcroft P, Huisman JS, Lehtinen S, Bouman JA, Althaus CL, Regoes RR, et al.
 COVID-19 infectivity profile correction. Swiss medical weekly. 2020. p. w20336.
*/
class InfecProfile_Ashcroft20 : public Base_InfecProfile
{
private:
	std::shared_ptr<Base_Distribution> gamma_distr_ptr;

	void defaults()
	{
		init("ashcroft20", {{"shape", 97.32}, {"scale", 0.268}, {"shift", 25.6}});
		gamma_distr_ptr = distribution_switchboard("gamma");
	}

public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Ashcroft20() : Base_InfecProfile() { defaults(); };
	InfecProfile_Ashcroft20(std::function<double(parID)> get) : Base_InfecProfile(get) { defaults(); }

	double pdf(const std::vector<double> &x) override
	{
		// load variables
		double t_OST = x[0] - x[1]; // t. from onset to transmission

		// update arguments
		gamma_distr_ptr->set_arg(0, get_arg(0));		   // shape
		gamma_distr_ptr->set_arg(1, get_arg(1));		   // scale
		gamma_distr_ptr->set_range({{"min", get_arg(2)}}); // UPDATE <===|
		double shift = get_arg(2);						   // shift argument

		return gamma_distr_ptr->pdf(shift + t_OST);
	}

	double cdf(const std::vector<double> &x) override
	{
		// load variables
		double t_OST = x[0] - x[1]; // t. from onset to transmission

		// update arguments
		gamma_distr_ptr->set_arg(0, get_arg(0));		   // shape
		gamma_distr_ptr->set_arg(1, get_arg(1));		   // scale
		gamma_distr_ptr->set_range({{"min", get_arg(2)}}); // UPDATE <===|
		double shift = get_arg(2);						   // shift argument

		return gamma_distr_ptr->cdf(shift + t_OST);
	}
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  MCMC : methods (user defined)  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the infectivity profile. */
std::unique_ptr<Base_InfecProfile> infecprofile_switchboard(const std::string &name,
															std::function<double(parID)> get)
{
	const str_vec options = {"ferretti20", "ashcroft20"};

	if (!tools::is_na(name))
	{
		if (name == "ferretti20")
			return std::make_unique<InfecProfile_Ferretti20>(std::move(get));

		else if (name == "ashcroft20")
			return std::make_unique<InfecProfile_Ashcroft20>(std::move(get));

		else
			throw std::invalid_argument(
				"Unknown infectivity profile '" + name + "'." + tools::find_closest(name, options));
	}

	return std::make_unique<Base_InfecProfile>();
}