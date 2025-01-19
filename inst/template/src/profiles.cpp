/*
 * profiles.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "distribution.h"
#include "mcmc_build.h"

// UPDATE: this should be merge into functions finally.
// It could be done by changing Arguments.

/* Base Class for Infectivity Profiles.
	NOTE : this classes access to parameters passing through MCMC class.
	This way I can check which parameters changes augmented dependencies. */
class MCMC::Impl::InfecProfile
{
public:
	//===== CONSTRUCTOR =====//
	InfecProfile(std::function<double(parID)> get = nullptr);
	virtual ~InfecProfile();

	//===== GETTER =====//
	double get_arg(ushort i);

	virtual double pdf(const dbl_vec &x);
	virtual double cdf(const dbl_vec &x);

protected:
	//===== SETTER =====//
	void _init_defaults_(const std::string &name,
						 const str_vec &arg_name,
						 const dbl_vec &arg_value);

private:
	std::string name = NULL_STRING;
	str_vec _arg_names_ = {};
	dbl_vec _values_ = {};
	std::vector<parID> _indexes_ = {};

	//----- Functional methods -----//
	double (InfecProfile::*_func_get_arg_)(ushort) = &InfecProfile::_internal_get_arg_;
	std::function<double(parID)> _external_param_getter_ = nullptr;

	double _internal_get_arg_(ushort i);
	double _external_get_arg_(ushort i);
};

MCMC::Impl::InfecProfile::InfecProfile(std::function<double(parID)> get) : _external_param_getter_(std::move(get)) {}
MCMC::Impl::InfecProfile::~InfecProfile() = default;

void MCMC::Impl::InfecProfile::_init_defaults_(const std::string &name,
											   const str_vec &arg_name,
											   const dbl_vec &arg_value)
{
	this->name = name;

	for (ushort c = 0; c < arg_name.size(); ++c)
	{
		this->_arg_names_.emplace_back(arg_name[c]);
		this->_values_.emplace_back(arg_value[c]);
	}
}

double MCMC::Impl::InfecProfile::pdf(const dbl_vec &x)
{
	UNUSED(x);
	throw std::invalid_argument(
		"Call to undefined method \"pdf\" of infectivity profile \"" + name + "\".");
}
double MCMC::Impl::InfecProfile::cdf(const dbl_vec &x)
{
	UNUSED(x);
	throw std::invalid_argument(
		"Call to undefined method \"cdf\" of infectivity profile \"" + name + "\".");
}

/* Return value indexed as 'i'. */
double MCMC::Impl::InfecProfile::get_arg(ushort i)
{
	return (this->*_func_get_arg_)(i);
}

/* Return internal value indexed as 'i'. */
double MCMC::Impl::InfecProfile::_internal_get_arg_(ushort i)
{
	return _values_[i];
}

/* Return external value indexed as 'i' with functional . */
double MCMC::Impl::InfecProfile::_external_get_arg_(ushort i)
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
class MCMC::Impl::InfecProfile_Ferretti20 : public MCMC::Impl::InfecProfile
{
public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Ferretti20(std::function<double(parID)> get = nullptr) : MCMC::Impl::InfecProfile(get)
	{
		_init_defaults_("ferretti20",
						{"alpha", "mu", "sigma", "tau"},
						{5.85, -4.0, 1.85, 5.8});
	}
	~InfecProfile_Ferretti20() = default;

	//*  Probability distribution function. *//
	double pdf(const dbl_vec &x) override
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
	double cdf(const dbl_vec &x) override
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

class InfecProfile_Ashcroft20 : public Base_InfecProfile
{
private:
	std::shared_ptr<Base_Distribution> gamma_distr_ptr;

public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Ashcroft20(std::function<double(parID)> get = nullptr) : Base_InfecProfile(get)
	{
		_init_defaults_("ashcroft20", {"shape", "scale", "shift"}, {97.32, 0.268, 25.6});
		gamma_distr_ptr = distribution_switchboard("gamma");
	}

	double pdf(const dbl_vec &x) override
	{
		// load variables
		double t_OST = x[0] - x[1]; // t. from onset to transmission

		// update arguments
		// UPDATE <===| these should just be assigned as pointers

		gamma_distr_ptr->set_arg(0, get_arg(0));				  // shape
		gamma_distr_ptr->set_arg(1, get_arg(1));				  // scale
		gamma_distr_ptr->set_range(get_arg(2), NULL_DOUBLE); // shift

		double shift = get_arg(2);								  // shift argument

		return gamma_distr_ptr->pdf(shift + t_OST);
	}

	double cdf(const dbl_vec &x) override
	{
		// load variables
		double t_OST = x[0] - x[1]; // t. from onset to transmission

		// update arguments
		// UPDATE <===| these should just be assigned as pointers

		gamma_distr_ptr->set_arg(0, get_arg(0));				  // shape
		gamma_distr_ptr->set_arg(1, get_arg(1));				  // scale
		gamma_distr_ptr->set_range(get_arg(2), NULL_DOUBLE); // shift
		double shift = get_arg(2);								  // shift argument

		return gamma_distr_ptr->cdf(shift + t_OST);
	}
};
*/

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////  MCMC : Impl : constructor  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

MCMC::Impl::Impl() : _uniform_distr_(0., 1.),
					 _infec_profile_(std::make_unique<MCMC::Impl::InfecProfile>())
{
#ifdef _OPENMP
	_max_n_available_threads = omp_get_max_threads();
#endif
}
MCMC::Impl::~Impl() = default;

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : infectivity profile  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the infectivity profile. */
void MCMC::Impl::assign_infectivity_profile(const std::string &name)
{
	auto lambda_getter = [this](const parID &p)
	{ return get_param(p); };

	// FLAG
	//	const str_vec options = {"ferretti20", "ashcroft20"};
	const str_vec options = {"ferretti20"};

	if (!tools::is_na(name))
	{
		if (name == "ferretti20")
			_infec_profile_ = std::make_unique<InfecProfile_Ferretti20>(std::move(lambda_getter));

		// else if (name == "ashcroft20") // FLAG
		//	return std::make_unique<InfecProfile_Ashcroft20>(std::move(lambda_getter));

		else
			throw std::invalid_argument(
				"Unknown infectivity profile \"" + name + "\"." + tools::find_closest(name, options));
	}

	_infec_profile_ = std::make_unique<InfecProfile>();
}

void MCMC::assign_infectivity_profile(const std::string &name)
{
	pBuild->assign_infectivity_profile(name);
}

/* Value of the instantaneous infectivity profile. */
double MCMC::infec_profile_inst(const dbl_vec &x) { return pBuild->_infec_profile_->pdf(x); }

/* Value of the cumulative infectivity profile. */
double MCMC::infec_profile_cumul(const dbl_vec &x) { return pBuild->_infec_profile_->cdf(x); }