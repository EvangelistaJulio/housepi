/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * mcmc_infec_profiles.cpp
 *
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
					  0., HOUSEPI_NA_DOUBLE, false, false);
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
/////////////////////////////////// ViralLoad_Cortier25 : class ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*

 */
class InfecProfile_Cortier25 : public Distribution
{
public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Cortier25()
	{
		init_defaults("cortier25",
					  {"gamma", "baseline_vl", "reference_vl"},
					  {2.5, -2, 0.5},
					  0., HOUSEPI_NA_DOUBLE, false, false);
	}
	~InfecProfile_Cortier25() = default;

	//===== SETTER =====//
	void init_reshaping(bool discrete) override final
	{
		init_min_zero_reshaping(discrete, false);
	}

protected:
	double raw_pdf(double x, const dbl_vec &opts) const override final
	{
		const double v = opts.at(0);  // (maximum VL - baseline VL)
		const double Tp = opts.at(1); // proliferation time
		const double Tc = opts.at(2); // clearence time

		const double g = get_arg(0);	// 'gamma' : shape
		const double vl0 = get_arg(1);	// 'gamma' : shape
		const double vref = get_arg(2); // 'gamma' : shape

		const double dt = x / Tp; // normalized time interval
		const double r = Tc / Tp; // ratio clear / proliferation time

		const double thres = std::max(0., -vl0 / v);

		if (thres >= 1)
			return 0.;

		else if (dt <= thres)
			return 0.;

		else if (dt <= 1.)
			return pow((vl0 + v * dt) / vref, g);

		else if (dt < 1. + r * (1. - thres))
			return pow((vl0 + v * (1. + r - dt) / r) / vref, g);

		else
			return 0.;
	}

	double raw_cdf(double x, const dbl_vec &opts) const override final
	{
		const double v = opts.at(0);  // (maximum VL - baseline VL)
		const double Tp = opts.at(1); // proliferation time
		const double Tc = opts.at(2); // clearence time

		const double g = get_arg(0);	// 'gamma' : shape
		const double vl0 = get_arg(1);	// 'gamma' : shape
		const double vref = get_arg(2); // 'gamma' : shape

		const double dt = x / Tp; // normalized time interval
		const double r = Tc / Tp; // ratio clear / proliferation time

		const double thres = std::max(0., -vl0 / v);

		const double den = v * (g + 1.) * pow(vref, g) / Tp;

		if (thres >= 1)
			return 0.;

		else if (dt <= thres)
			return 0.;

		else if (dt <= 1.)
			return (pow(vl0 + v * dt, g + 1.) - pow(vl0 + v * thres, g + 1.)) / den;

		else if (dt < 1. + r * (1. - thres))
		{
			const double F1 = pow(vl0 + v, g + 1.) - pow(vl0 + v * thres, g + 1.);
			return (F1 + r * (pow(vl0 + v, g + 1) - pow(vl0 + v * (1. + r - dt) / r, g + 1))) / den;
		}
		else
		{
			const double F1 = pow(vl0 + v, g + 1.) - pow(vl0 + v * thres, g + 1.);
			return F1 * (1 + r) / den;
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
////////////////////////////////// InfecProfile_Ashcroft20 : class ////////////////////////////////
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
					  0., HOUSEPI_NA_DOUBLE, false, false);

		// wrapper around gamma distribution
		gamma_distr_ptr = distribution_switch("gamma");
		gamma_distr_ptr->assign_attr_pointer("shape", ATTR_ARGUMENT, get_ptr_to_arg(0));
		gamma_distr_ptr->assign_attr_pointer("scale", ATTR_ARGUMENT, get_ptr_to_arg(1));
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
void MCMC::Impl::InfecProfile::_assign_profile(const std::string &name)
{
	const str_vec options = {"ferretti20", "ashcroft20", "cortier25"};

	if (name == "ferretti20")
		profile = std::make_shared<InfecProfile_Ferretti20>();

	else if (name == "ashcroft20")
		profile = std::make_shared<InfecProfile_Ashcroft20>();

	else if (name == "cortier25")
		profile = std::make_shared<InfecProfile_Cortier25>();

	else
	{
		try
		{
			profile = distribution_switch(name);
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

/* Adjust lists of parameters to be inferred making parameter as dependent on cache. */
void MCMC::Impl::make_param_cache_dependent(const std::string &par_name)
{
	// remove a parameter from standard if it was defined as such
	auto it = std::find(param_standard_infer_list.begin(), param_standard_infer_list.end(), par_name);
	if (it != param_standard_infer_list.end())
		param_standard_infer_list.erase(it);

	// add the parameter as cache
	it = std::find(param_cache_infer_list.begin(), param_cache_infer_list.end(), par_name);
	if (it == param_cache_infer_list.end())
		param_cache_infer_list.push_back(par_name);
}

/* Assign arguments to infectivity profile. */
void MCMC::Impl::InfecProfile::_assign_args(const str_vec &keys,
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
			if (pImpl.param_map.at(value).is_constant())
				profile->assign_attr_scalar(key, ATTR_ARGUMENT, pImpl.get_param(value));
			else
			{
				profile->assign_attr_pointer(key, ATTR_ARGUMENT, pImpl.get_param_pointer(value));
				pImpl.make_param_cache_dependent(value);
			}
		}
		else if (type == TYPE_SCALAR)
			profile->assign_attr_scalar(key, ATTR_ARGUMENT, tools::stod(value, "argument with scalar type"));

		else
			throw std::invalid_argument(
				"Unsopported type \"" + type + "\" found while loading \"" + value + "\" for pointer value.");
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : InfecProfile : assign range  ////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Assign maximum of infectivity profile. */
void MCMC::Impl::InfecProfile::_assign_range(const std::string &range_name,
											 const std::string &max_value,
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
		if (!pImpl.param_map.at(max_value).is_constant())
		{
			profile->assign_attr_pointer(range_name, ATTR_RANGE, pImpl.get_param_pointer(max_value));
			pImpl.make_param_cache_dependent(max_value);
		}
		else
			profile->assign_attr_scalar(range_name, ATTR_RANGE, pImpl.get_param(max_value));
	}

	else if (max_at == TYPE_SCALAR)
		profile->assign_attr_scalar(range_name, ATTR_RANGE, tools::stod(max_value, "argument with scalar type"));

	else
		throw std::invalid_argument(
			"Unsopported type \"" + max_at + "\" found while loading \"" + max_value + "\" for pointer value.");
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////  MCMC : Impl : InfecProfile : assign template  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Function to assign a template to infectivity profile. */
void MCMC::Impl::InfecProfile::assign_template(std::map<std::string, std::string> &x)
{
	if (deflist.empty())
	{
		if (!x.count("name"))
			throw std::invalid_argument("Missing name for infectivity profile template.");

		for (const std::string &key : _headers)
			if (!x.count(key))
				x.emplace(key, HOUSEPI_NA_STRING);

		if (tools::is_na(x.at("discr")))
			x["discr"] = "false";

		deflist = x;
		build();
	}
}

/* Store values of the infectivity profile. */
void MCMC::Impl::InfecProfile::init_lookup()
{
	if (profile->is_discrete())
	{
		bool is_suitable = true;
		// UPDATE : check if requires optional arguments
		try
		{
			UNUSED(__default_inst(profile->get_min() + 1, {}));
		}
		catch (std::out_of_range &e)
		{
			is_suitable = false;
		}

		if (is_suitable)
		{
			if (tools::is_na(profile->get_max()))
				_size_lut = DEFAULT_MAX_SIZE;
			else
				_size_lut = profile->get_max() + 1;

			_lookup_inst.resize(_size_lut);
			_lookup_cumul.resize(_size_lut);
			_lookup_buffer_inst.resize(_size_lut);
			_lookup_buffer_cumul.resize(_size_lut);

			update_lookup();

			_do_inst = &InfecProfile::__lookup_inst;
			_do_cumul = &InfecProfile::__lookup_cumul;
		}
	}
}

/* Restore buffered lookup table. */
void MCMC::Impl::InfecProfile::restore_lookup()
{
	_lookup_inst.swap(_lookup_buffer_inst);
	_lookup_cumul.swap(_lookup_buffer_cumul);
}

/* Store values of the infectivity profile. */
void MCMC::Impl::InfecProfile::update_lookup()
{
	restore_lookup();

	// Set infectivity to 0 on the same day
	_lookup_inst[0] = 0.;
	_lookup_cumul[0] = 0.;

	for (size_t i = 1; i < _size_lut; ++i)
	{
		double x = __default_inst(i, {});
		_lookup_inst[i] = x;
		_lookup_cumul[i] = x + _lookup_cumul[i - 1];
	}
}

/* Build infectivity profile distribution . */
void MCMC::Impl::InfecProfile::build()
{
	try
	{
		if (!deflist.empty())
		{
			// check entries
			tools::check_for_unused_entries(tools::headers_from_deflist(deflist),
											_headers);
			bool is_discrete = true;

			//----- Name -----//
			if ((deflist.count("name")) && (!tools::is_na(deflist.at("name"))))
				_assign_profile(deflist.at("name"));
			else
				throw std::invalid_argument("missing name.");

			//----- Discrete -----//
			if (!profile->is_discrete()) // defined continuos
			{
				if (deflist.count("discr"))
					is_discrete = tools::stobool(deflist.at("discr"));
				else
					throw std::invalid_argument("not a valid discrete.");
			}

			//----- Arguments -----//
			if (deflist.count("arg_key") && deflist.count("arg_val") && deflist.count("arg_at"))
				_assign_args(tools::parse_list_args(deflist.at("arg_key")),
							 tools::parse_list_args(deflist.at("arg_val")),
							 tools::parse_list_args(deflist.at("arg_at")));
			else
				throw std::invalid_argument("not a valid list of arguments.");

			//----- Range minimum -----//
			if (deflist.count("min") && deflist.count("min_at"))
			{
				if (tools::is_na(deflist.at("min")) && tools::is_na(deflist.at("min_at")))
					_assign_range("min", (is_discrete) ? "1" : "0", TYPE_SCALAR);
				else
					_assign_range("min", deflist.at("min"), deflist.at("min_at"));

				// check infectivity profile minimum is correct
				if (is_discrete)
				{
					if (profile->get_min() < 1.)
						throw std::invalid_argument(
							"Minimum of discrete infectivity profile must be >= 1.");
				}
				else if (profile->get_min() < 0.)
					throw std::invalid_argument(
						"Minimum of continuos infectivity profile must be >= 0.");
			}
			else
				throw std::invalid_argument("not a valid minimum.");

			//----- Range maximum -----//
			if (deflist.count("max") && deflist.count("max_at"))
				_assign_range("max", deflist.at("max"), deflist.at("max_at"));
			else
				throw std::invalid_argument("not a valid maximum.");

			//----- Reshaping -----//
			if (deflist.count("discr"))
				profile->init_reshaping(tools::stobool(deflist.at("discr")));
			else
				throw std::invalid_argument("not a valid discrete.");

			// if discrete and suitable initializes the lookup table
			init_lookup();
		}
	}
	catch (const std::exception &e)
	{
		tools::error({"While building infectivity profile definition:",
					  std::string(e.what())});
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  MCMC : Impl : InfecProfile : getter methods  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Check if a lookup table exists. */
bool MCMC::Impl::InfecProfile::is_lookup() const { return !_lookup_inst.empty(); }

/* Default method for instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::__default_inst(double t, const dbl_vec &opts) const
{
	return profile->pdf(t, opts);
}
/* Default method for cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::__default_cumul(double t, const dbl_vec &opts) const
{
	return profile->cdf(t, opts);
}

/* Lookup method for discrete instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::__lookup_inst(double t, const dbl_vec &opts) const
{
	UNUSED(opts);
	if ((t > 0) && (t < _size_lut))
		return _lookup_inst.at(static_cast<size_t>(t));
	else
		return 0;
}
/* Lookup method for discrete cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::__lookup_cumul(double t, const dbl_vec &opts) const
{
	UNUSED(opts);
	if (t < 0)
		return 0.;
	else if (t >= _size_lut)
		return 1.;
	else
		return _lookup_cumul.at(static_cast<size_t>(t));
}

/* Value of the instantaneous infectivity profile. */
double MCMC::Impl::InfecProfile::inst(double t, const dbl_vec &opts) const
{
	if (profile->valid(t))
		return (this->*_do_inst)(t, opts);
	else
		return 0.;
}
/* Value of the cumulative infectivity profile. */
double MCMC::Impl::InfecProfile::cumul(double t, const dbl_vec &opts) const
{
	if (profile->below(t))
		return 0.;
	else if (profile->above(t))
		return 1.;
	else
		return (this->*_do_cumul)(t, opts);
}

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  MCMC : Impl : update\restore lookup tables  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* If infectivity profile is a lookup table, store its values. */
void MCMC::Impl::update_infec_profile()
{
	if (_infec_profile.is_lookup())
		_infec_profile.update_lookup();
}

/* If infectivity profile is a lookup table, restore its buffered values. */
void MCMC::Impl::restore_infec_profile()
{
	if (_infec_profile.is_lookup())
		_infec_profile.restore_lookup();
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////  MCMC : infectivity profile  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Virtual method for user initialization of infectivity profile. */
void MCMC::init_infec_profile() {};

/* Assign default infectivity profile during user model initialization. */
void MCMC::assign_template_infec_profile(std::map<std::string, std::string> x)
{
	pImpl->_infec_profile.assign_template(x);
}

/* Value of the instantaneous infectivity profile. */
double MCMC::infec_profile_inst(double t, const dbl_vec &opts) const
{
	return pImpl->_infec_profile.inst(t, opts);
}

/* Value of the cumulative infectivity profile. */
double MCMC::infec_profile_cumul(double t, const dbl_vec &opts) const
{
	return pImpl->_infec_profile.cumul(t, opts);
}
