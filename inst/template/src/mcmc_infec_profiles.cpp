/*
 * mcmc_infec_profiles.cpp
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 */

#include "distribution.h"
#include "mcmc_build.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

///////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// InfecProfile_Ferretti20 : class ////////////////////////////////
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
class InfecProfile_Ferretti20 : public Distribution
{
public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Ferretti20()
	{
		_init_defaults_("ferretti20",
						{"shape", "loc", "scale", "incub"},
						{5.85, -4.0, 1.85, 5.8},
						NULL_DOUBLE, NULL_DOUBLE, false);
	}
	~InfecProfile_Ferretti20() = default;

protected:
	double _pdf_(double x, const dbl_vec &opts) const override
	{
		double i = opts.at(0); // incubation period
		double a = get_arg(0); // 'alpha' : shape
		double m = get_arg(1); // 'mu' : location
		double s = get_arg(2); // 'sigma' : scale
		double T = get_arg(3); // 'tau' : mean incubation period
		double tOST = x - i;   // t. from onset to transmission

		double norm_pdf = (s / a) * (1. + (i / T - 1.) * pow(1 + exp(m / s), -a) - (i / T) * pow(1. + exp((m + T) / s), -a));

		if (tOST < 0)
			tOST *= (T / i);

		return exp((m - tOST) / s) * pow(1. + exp((m - tOST) / s), -a - 1.) / norm_pdf;
	}

	double _cdf_(double x, const dbl_vec &opts) const override
	{
		double i = opts.at(0); // incubation period
		double a = get_arg(0); // 'alpha' : shape
		double m = get_arg(1); // 'mu' : location
		double s = get_arg(2); // 'sigma' : scale
		double T = get_arg(3); // 'tau' : mean incubation period
		double tOST = x - i;   // t. from onset to transmission

		double norm_cdf = (1. + (i / T - 1.) * pow(1 + exp(m / s), -a) - (i / T) * (pow(1. + exp((m + T) / s), -a)));

		if (tOST < 0) // Integrate from -incub_period to tOST
			return (i / T) * (pow(1. + exp((m - tOST * T / i) / s), -a) - pow(1. + exp((m + T) / s), -a)) / norm_cdf;

		else
		{
			// Integrate from -incub to 0
			double tmp = (i / T) * (pow(1. + exp(m / s), -a) - pow(1. + exp((m + T) / s), -a));
			// Integrate from 0 to t_ost
			return (tmp + pow(1. + exp((m - tOST) / s), -a) - pow(1. + exp(m / s), -a)) / norm_cdf;
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
class InfecProfile_Ashcroft20 : public Distribution
{
private:
	std::shared_ptr<Distribution> gamma_distr_ptr;

public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Ashcroft20()
	{
		_init_defaults_("ashcroft20",
						{"shape", "scale", "shift"},
						{97.32, 0.268, 25.6},
						NULL_DOUBLE, NULL_DOUBLE, false);

		gamma_distr_ptr = distribution_switchboard("gamma");
		gamma_distr_ptr->assign_pointer_arg("shape", get_ptr_to_arg(0));
		gamma_distr_ptr->assign_pointer_arg("scale", get_ptr_to_arg(1));
		gamma_distr_ptr->assign_pointer_range_min(get_ptr_to_arg(2));
	}
	~InfecProfile_Ashcroft20() = default;

protected:
	double _pdf_(double x, const dbl_vec &opts) const override
	{
		double i = opts.at(0); // incubation period
		double tOST = x - i;   // t. from onset to transmission
		return gamma_distr_ptr->pdf(tOST);
	}

	double _cdf_(double x, const dbl_vec &opts) const override
	{
		double i = opts.at(0); // incubation period
		double tOST = x - i;   // t. from onset to transmission
		return gamma_distr_ptr->cdf(tOST);
	}
};

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  MCMC : Impl : infectivity profile  ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign the infectivity profile. */
void MCMC::Impl::switchboard_infec_profile(const std::string &name)
{
	const str_vec options = {"ferretti20", "ashcroft20"};

	if (name == "ferretti20")
		_infec_profile_ = std::make_shared<InfecProfile_Ferretti20>();

	else if (name == "ashcroft20")
		_infec_profile_ = std::make_shared<InfecProfile_Ashcroft20>();

	else
	{
		try
		{
			_infec_profile_ = distribution_switchboard(name);
		}
		catch (std::exception &e)
		{
			throw std::invalid_argument(
				e.what() + tools::find_closest(name, options));
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC : Impl : assign infectivity profile arguments  ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign arguments to infectivity profile. */
void MCMC::Impl::assign_infec_profile_args(const str_vec &keys,
										   const str_vec &vals,
										   const str_vec &typs)
{
	// check size validity of keys, values and type
	if ((keys.size() != vals.size()) || (typs.size() != vals.size()))
		throw std::invalid_argument(
			"Mismatching definitions in arguments.");

	// loop over all arguments
	for (size_t x = 0; x < keys.size(); ++x)
	{
		std::string key = keys[x];
		std::string value = vals[x];
		std::string type = typs[x];

		if (tools::is_na(value) && tools::is_na(type))
			; // pass

		else if (!tools::is_na(value) && tools::is_na(type))
			throw std::invalid_argument(
				"Impossible to interpret default value of  \"" + value + "\". Please provide a type.");

		else if (tools::is_na(value) && !tools::is_na(type))
			tools::warning(
				"Found a type without a paired value. Falling back to default.");

		else if (type == TYPE_PARAM)
		{
			if (param_map.at(value).is_constant())
				_infec_profile_->assign_scalar_arg(key, get_param(value));
			else
				_infec_profile_->assign_pointer_arg(key, get_param_pointer(value));
		}
		else if (type == TYPE_SCALAR)
			_infec_profile_->assign_scalar_arg(key, tools::stod(value, "argument with scalar type"));

		else
			throw std::invalid_argument(
				"Unsopported type \"" + type + "\" found while loading \"" + value + "\" for pointer value.");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////  MCMC : Impl : assign infectivity profile range  //////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Convert a std::string with type-dictionary to std::string-to-double */
void MCMC::Impl::assign_infec_profile_range(const std::string &min_value,
											const std::string &min_type,
											const std::string &max_value,
											const std::string &max_type)
{
	//-------------------//
	//----- Minimum -----//
	//-------------------//

	if (tools::is_na(min_value) && tools::is_na(min_type))
		; // pass

	else if (!tools::is_na(min_value) && tools::is_na(min_type))
		throw std::invalid_argument(
			"Impossible to interpret default value of  \"" + min_value + "\". Please provide a type.");

	else if (tools::is_na(min_value) && !tools::is_na(min_type))
		tools::warning(
			"Found a type without a paired value. Falling back to default.");

	else if (min_value == TYPE_PARAM)
	{
		if (param_map.at(min_value).is_constant())
			_infec_profile_->assign_scalar_range_min(get_param(min_value));
		else
			_infec_profile_->assign_pointer_range_min(get_param_pointer(min_value));
	}

	else if (min_type == TYPE_SCALAR)
		_infec_profile_->assign_scalar_range_min(tools::stod(min_value, "scalar type argument"));

	else
		throw std::invalid_argument(
			"Unsopported type \"" + min_type + "\" found while loading \"" + min_value + "\" for pointer value.");

	//-------------------//
	//----- Maximum -----//
	//-------------------//

	if (tools::is_na(max_value) && tools::is_na(max_type))
		; // pass

	else if (!tools::is_na(max_value) && tools::is_na(max_type))
		throw std::invalid_argument(
			"Impossible to interpret default value of  \"" + max_value + "\". Please provide a type.");

	else if (tools::is_na(max_value) && !tools::is_na(max_type))
		tools::warning(
			"Found a type without a paired value. Falling back to default.");

	else if (max_type == TYPE_PARAM)
	{
		if (param_map.at(max_value).is_constant())
			_infec_profile_->assign_scalar_range_max(get_param(max_value));
		else
			_infec_profile_->assign_pointer_range_max(get_param_pointer(max_value));
	}

	else if (max_type == TYPE_SCALAR)
		_infec_profile_->assign_scalar_range_max(tools::stod(max_value, "argument with scalar type"));

	else
		throw std::invalid_argument(
			"Unsopported type \"" + max_type + "\" found while loading \"" + max_value + "\" for pointer value.");
}


/* Read Setup file and produce a definition list. */
void MCMC::Impl::read_infec_profile(const std::string &pathtofile)
{
    tools::read_csv_list(pathtofile, _infec_profile_deflist_);
    // note: build_infec_profile is optionally called by the user
    _report_pathtofile_["infec_profile"] = pathtofile;
}

/* Build infectivity profile distribution . */
void MCMC::Impl::build_infec_profile()
{
	try
	{
		// reference to infectivity profile definition list
		std::map<std::string, std::string> &deflist = _infec_profile_deflist_;

		if (!deflist.empty())
		{
			//===== Check that all entries are correct =====//
			const str_vec options = {"name", "discrete",
									 "arg_name", "arg_value", "arg_type",
									 "min", "min_type",
									 "max", "max_type"};
			str_vec entries = {};
			for (const auto &pair : deflist)
				entries.push_back(pair.first);
			tools::check_unused_entries(entries, options);

			//=====================================================//
			//===== ASSIGN INFECTIVITY PROFILE CONFIGURATIONS =====//
			//=====================================================//

			//----- Name -----//
			if ((deflist.count("name")) && (!tools::is_na(deflist.at("name"))))
				switchboard_infec_profile(deflist.at("name"));
			else
				throw std::invalid_argument("Provide a valid name.");

			//----- Arguments -----//
			if ((deflist.count("arg_name")) && (deflist.count("arg_value")) && (deflist.count("arg_type")))
				assign_infec_profile_args(tools::parse_list_args(deflist.at("arg_name")),
										  tools::parse_list_args(deflist.at("arg_value")),
										  tools::parse_list_args(deflist.at("arg_type")));
			else
				throw std::invalid_argument("Provide a valid list of arguments.");

			//----- Range -----//
			if ((deflist.count("min")) && (deflist.count("min_type")))
				if ((deflist.count("max")) && (deflist.count("max_type")))
					assign_infec_profile_range(deflist.at("min"), deflist.at("min_type"),
											   deflist.at("max"), deflist.at("max_type"));
				else
					throw std::invalid_argument("Provide a valid maximum.");
			else
				throw std::invalid_argument("Provide a valid minimum.");

			//----- Reshaping -----//
			if (deflist.count("discrete"))
				_infec_profile_->init_reshaping(tools::stobool(deflist.at("discrete")));
			else
				throw std::invalid_argument("Provide a valid minimum.");
		}
	}
	catch (const std::exception &e)
	{
		tools::error(
			"While building infectivity profile definition: " + std::string(e.what()));
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : infectivity profile  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Call to infectivity profile during user model initialization. */
void MCMC::assign_infec_profile(const std::string &default_name)
{
	if (pImpl->_infec_profile_deflist_.empty())
		// default case
		pImpl->switchboard_infec_profile(default_name);
	else
		// user defined case
		pImpl->build_infec_profile();
}

/* Value of the instantaneous infectivity profile. */
double MCMC::infec_profile_inst(double gener_time, const dbl_vec &opts) const
{
	return pImpl->_infec_profile_->pdf(gener_time, opts);
}

/* Value of the cumulative infectivity profile. */
double MCMC::infec_profile_cumul(double gener_time, const dbl_vec &opts) const
{
	return pImpl->_infec_profile_->cdf(gener_time, opts);
}
