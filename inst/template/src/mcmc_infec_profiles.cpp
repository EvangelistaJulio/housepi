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
		init_defaults("ferretti20",
					  {"shape", "loc", "scale", "incub"},
					  {5.85, -4.0, 1.85, 5.8},
					  0., NULL_DOUBLE, false, false);
	}
	~InfecProfile_Ferretti20() = default;

	//===== SETTER =====//
	void init_reshaping(bool discrete) override final
	{
		init_min_zero_reshaping(discrete, false);
	}

protected:
	double raw_pdf(double x, const dbl_vec &opts) const override final
	{
		double i = opts.at(0); // incubation period
		double a = get_arg(0); // 'alpha' : shape
		double m = get_arg(1); // 'mu' : location
		double s = get_arg(2); // 'sigma' : scale
		double T = get_arg(3); // 'tau' : mean incubation period
		double tOST = x - i;   // time from onset to transmission

		double norm_pdf = (s / a) * (1. + (i / T - 1.) * pow(1 + exp(m / s), -a) - (i / T) * pow(1. + exp((m + T) / s), -a));

		if (tOST < 0)
			tOST *= (T / i);

		return exp((m - tOST) / s) * pow(1. + exp((m - tOST) / s), -a - 1.) / norm_pdf;
	}

	double raw_cdf(double x, const dbl_vec &opts) const override final
	{
		double i = opts.at(0); // incubation period
		double a = get_arg(0); // 'alpha' : shape
		double m = get_arg(1); // 'mu' : location
		double s = get_arg(2); // 'sigma' : scale
		double T = get_arg(3); // 'tau' : mean incubation period
		double tOST = x - i;   // time from onset to transmission

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
public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Ashcroft20()
	{
		init_defaults("ashcroft20",
					  {"shape", "scale", "shift"},
					  {97.32, 0.268, 25.6},
					  0., NULL_DOUBLE, false, false);

		// wrapper around gamma distribution
		gamma_distr_ptr = distribution_switchboard("gamma");
		gamma_distr_ptr->assign_pointer_arg("shape", get_ptr_to_arg(0));
		gamma_distr_ptr->assign_pointer_arg("scale", get_ptr_to_arg(1));
	}
	~InfecProfile_Ashcroft20() = default;

	//===== SETTER =====//
	void init_reshaping(bool discrete) override final
	{
		init_min_zero_reshaping(discrete, false);
	}

protected:
	std::shared_ptr<Distribution> gamma_distr_ptr;

	double raw_pdf(double x, const dbl_vec &opts) const override final
	{
		double i1 = opts.at(0);		 // incubation period infector
		double i2 = opts.at(1);		 // incubation period infectee
		double serial = x + i2 - i1; // serial time
		return gamma_distr_ptr->pdf(serial + get_arg(2));
	}

	double raw_cdf(double x, const dbl_vec &opts) const override final
	{
		double i1 = opts.at(0);		 // incubation period infector
		double i2 = opts.at(1);		 // incubation period infectee
		double serial = x + i2 - i1; // serial time
		return gamma_distr_ptr->cdf(serial + get_arg(2));
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
	const str_vec options = {"ferretti20"};

	if (name == "ferretti20")
		_infec_profile_ = std::make_shared<InfecProfile_Ferretti20>();

	else
	{
		try
		{
			_infec_profile_ = distribution_switchboard(name);
		}
		catch (std::exception &e)
		{
			std::string prev = e.what();
			if (prev.back() == '?') // found closest amongst distributions
				throw std::invalid_argument(prev);

			else // search between infectivity profiles
			{
				std::string curr = tools::find_closest(name, options, 3);
				if (curr.back() == '?') // found closest amongst infec profiles
					throw std::invalid_argument(
						"Unknown infectivity profile \"" + name + "\"." + curr);

				else // combined error messages
					throw std::invalid_argument(
						prev + " Implemented infectivity profile(s): " + tools::show_list(options) + ".");
			}
		}
	}
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////  MCMC : Impl : assign infectivity profile arguments  ////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Adjust lists of parameters to be inferred */
void MCMC::Impl::fix_param_to_be_inferred_list(const std::string &par_name)
{
	// remove a parameter from standard if it was defined as such
	auto it = std::find(param_standard_infer_list.begin(), param_standard_infer_list.end(), par_name);
	if (it != param_standard_infer_list.end())
		param_standard_infer_list.erase(it);

	// add the parameter as extra
	it = std::find(param_extra_infer_list.begin(), param_extra_infer_list.end(), par_name);
	if (it == param_extra_infer_list.end())
		param_extra_infer_list.push_back(par_name);
}

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
			if (!param_map.at(value).is_constant())
			{
				_infec_profile_->assign_pointer_arg(key, get_param_pointer(value));
				fix_param_to_be_inferred_list(value);
			}
			else
				_infec_profile_->assign_scalar_arg(key, get_param(value));
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

/* Assign minimum of infectivity profile. */
void MCMC::Impl::assign_infec_profile_min(bool is_discrete,
										  const std::string &min_value,
										  const std::string &min_at)
{
	if (tools::is_na(min_value) && tools::is_na(min_at))
		_infec_profile_->assign_scalar_range_min((is_discrete) ? 1 : 0);

	else if (!tools::is_na(min_value) && tools::is_na(min_at))
		throw std::invalid_argument(
			"Impossible to interpret type of value  \"" + min_value + "\". Please provide a type.");

	else if (tools::is_na(min_value) && !tools::is_na(min_at))
		tools::warning(
			"Found a type without a paired value. Falling back to default.");

	else if (min_value == TYPE_PARAM)
	{
		if (!param_map.at(min_value).is_constant())
		{
			_infec_profile_->assign_pointer_range_min(get_param_pointer(min_value));
			fix_param_to_be_inferred_list(min_value);
		}
		else
			_infec_profile_->assign_scalar_range_min(get_param(min_value));
	}

	else if (min_at == TYPE_SCALAR)
		_infec_profile_->assign_scalar_range_min(tools::stod(min_value, "scalar type argument"));

	else
		throw std::invalid_argument(
			"Unsopported type \"" + min_at + "\" found while loading \"" + min_value + "\" for pointer value.");

	// check infectivity profile minimum is correct
	if (is_discrete)
	{
		if (_infec_profile_->get_min() < 1.)
			throw std::invalid_argument(
				"Minimum of discrete infectivity profile must be >= 1.");
	}
	else if (_infec_profile_->get_min() < 0.)
		throw std::invalid_argument(
			"Minimum of continuos infectivity profile must be >= 0.");
}

/* Assign maximum of infectivity profile. */
void MCMC::Impl::assign_infec_profile_max(const std::string &max_value,
										  const std::string &max_at)
{
	if (tools::is_na(max_value) && tools::is_na(max_at))
		; // pass

	else if (!tools::is_na(max_value) && tools::is_na(max_at))
		throw std::invalid_argument(
			"Impossible to interpret type of value  \"" + max_value + "\". Please provide a type.");

	else if (tools::is_na(max_value) && !tools::is_na(max_at))
		tools::warning(
			"Found a type without a paired value. Falling back to default.");

	else if (max_at == TYPE_PARAM)
	{
		if (!param_map.at(max_value).is_constant())
		{
			_infec_profile_->assign_pointer_range_max(get_param_pointer(max_value));
			fix_param_to_be_inferred_list(max_value);
		}
		else
			_infec_profile_->assign_scalar_range_max(get_param(max_value));
	}

	else if (max_at == TYPE_SCALAR)
		_infec_profile_->assign_scalar_range_max(tools::stod(max_value, "argument with scalar type"));

	else
		throw std::invalid_argument(
			"Unsopported type \"" + max_at + "\" found while loading \"" + max_value + "\" for pointer value.");
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
			// check entries
			str_vec headers = tools::headers_from_deflist(deflist);
			tools::check_for_unused_entries(headers, infec_profile_file_headers);
			bool is_discrete = true;

			//----- Name -----//
			if ((deflist.count("name")) && (!tools::is_na(deflist.at("name"))))
				switchboard_infec_profile(deflist.at("name"));
			else
				throw std::invalid_argument("missing name.");

			//----- Discrete -----//
			if (!_infec_profile_->is_discrete()) // defined continuos
			{
				if (deflist.count("discrete"))
					is_discrete = tools::stobool(deflist.at("discrete"));
				else
					throw std::invalid_argument("not a valid discrete.");
			}

			//----- Arguments -----//
			if ((deflist.count("arg_name")) && (deflist.count("arg_value")) && (deflist.count("arg_at")))
				assign_infec_profile_args(tools::parse_list_args(deflist.at("arg_name")),
										  tools::parse_list_args(deflist.at("arg_value")),
										  tools::parse_list_args(deflist.at("arg_at")));
			else
				throw std::invalid_argument("not a valid list of arguments.");

			//----- Range minimum -----//
			if ((deflist.count("min")) && (deflist.count("min_at")))
				assign_infec_profile_min(is_discrete, deflist.at("min"), deflist.at("min_at"));
			else
				throw std::invalid_argument("not a valid minimum.");

			//----- Range maximum -----//
			if ((deflist.count("max")) && (deflist.count("max_at")))
				assign_infec_profile_max(deflist.at("max"), deflist.at("max_at"));
			else
				throw std::invalid_argument("not a valid maximum.");

			//----- Reshaping -----//
			if (deflist.count("discrete"))
				_infec_profile_->init_reshaping(tools::stobool(deflist.at("discrete")));
			else
				throw std::invalid_argument("not a valid discrete.");
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

/* Virtual method for user initialization of infectivity profile. */
void MCMC::init_infec_profile() {};

/* Assign default infectivity profile during user model initialization. */
void MCMC::assign_default_infec_profile(const std::string &default_name)
{
	if (pImpl->_infec_profile_deflist_.empty())
		pImpl->switchboard_infec_profile(default_name);
}

/* Value of the instantaneous infectivity profile. */
double MCMC::infec_profile_inst(double gener_time, const dbl_vec &opts) const
{
	if (gener_time < 0)
		std::overflow_error(
			"Negative generation time passed infectivity profile.");

	return pImpl->_infec_profile_->pdf(gener_time, opts);
}

/* Value of the cumulative infectivity profile. */
double MCMC::infec_profile_cumul(double gener_time, const dbl_vec &opts) const
{
	if (gener_time < 0)
		std::overflow_error(
			"Negative generation time passed infectivity profile.");

	return pImpl->_infec_profile_->cdf(gener_time, opts);
}
