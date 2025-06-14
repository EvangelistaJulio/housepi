/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * infec_profiles.cpp
 *
 */

#include "distribution.h"

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
  0) "shape";
  1) "loc";
  2) "scale";
  3) "mean_incub" (mean incubation period);
  4) "incub_period" (individual incubation period);
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
		: Distribution("ferretti20",
					   {"shape", "loc", "scale", "mean_incub", "incub_period"},
					   {5.85, -4.0, 1.85, 5.8, 5.8},
					   0., HOUSEPI_NA_DOUBLE, false, false) {}
	~InfecProfile_Ferretti20() = default;

	//===== SETTER =====//
	void init_reshaping(bool discrete) override final
	{
		init_min_zero_reshaping(discrete, false);
	}

protected:
	double raw_pdf(double x) const override final
	{
		double a = get_arg(0); // 'shape'
		double m = get_arg(1); // 'loc'
		double s = get_arg(2); // 'scale'
		double T = get_arg(3); // 'mean_incub'
		double i = get_arg(4); // 'incub_period'
		double tOST = x - i;   // time from onset to transmission

		double norm_pdf = (s / a) * (1. + (i / T - 1.) * pow(1 + exp(m / s), -a) - (i / T) * pow(1. + exp((m + T) / s), -a));

		if (tOST < 0)
			tOST *= (T / i);

		return exp((m - tOST) / s) * pow(1. + exp((m - tOST) / s), -a - 1.) / norm_pdf;
	}

	double raw_cdf(double x) const override final
	{
		double a = get_arg(0); // 'shape'
		double m = get_arg(1); // 'loc'
		double s = get_arg(2); // 'scale'
		double T = get_arg(3); // 'mean_incub'
		double i = get_arg(4); // 'incub_period'
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
/////////////////////////////////// InfecProfile_Cortier25 : class ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/*

 */
class InfecProfile_Cortier25 : public Distribution
{
public:
	//===== CONSTRUCTOR =====//
	InfecProfile_Cortier25()
		: Distribution("cortier25",
					   {"shape", "refer_vl", "base_vl", "max_vl",
						"prolif_period", "clear_period"},
					   {2.5, 0.5, -2, 8, 2.4, 3.1},
					   0., HOUSEPI_NA_DOUBLE, false, false) {}
	~InfecProfile_Cortier25() = default;

	//===== SETTER =====//
	void init_reshaping(bool discrete) override final
	{
		// FIXME: I don't want truncation renormalization since no distribution
		init_min_zero_reshaping(discrete, false);
	}

protected:
	double raw_pdf(double x) const override final
	{
		const double s = get_arg(0);				// 'shape'
		const double vref = get_arg(1);				// 'refer_vl'
		const double vl0 = get_arg(2) / vref;		// 'base_vl' / 'refer_vl'
		const double v = (get_arg(3) - vl0) / vref; // ('max_vl' - 'base_vl') / 'refer_vl'
		const double Tp = get_arg(4);				// 'prolif_period'
		const double Tc = get_arg(5);				// 'clear_period'

		const double dt = x / Tp; // normalized time interval
		const double r = Tc / Tp; // ratio clear / proliferation time

		const double thres = std::max(0., -vl0 / v); // threshold VL > 0

		if (thres >= 1)
			return 0.;

		else if (dt <= thres)
			return 0.;

		else if (dt <= 1.)
			return pow(vl0 + v * dt, s);

		else if (dt < 1. + r * (1. - thres))
			return pow(vl0 + v * (1. + r - dt) / r, s);

		else
			return 0.;
	}

	double raw_cdf(double x) const override final
	{
		const double s = get_arg(0);				// 'shape'
		const double vref = get_arg(1);				// 'refer_vl'
		const double vl0 = get_arg(2) / vref;		// 'base_vl' / 'refer_vl'
		const double v = (get_arg(3) - vl0) / vref; // ('max_vl' - 'base_vl') / 'refer_vl'
		const double Tp = get_arg(4);				// 'prolif_period'
		const double Tc = get_arg(5);				// 'clear_period'

		const double dt = x / Tp; // normalized time interval
		const double r = Tc / Tp; // ratio clear / proliferation time

		const double thres = std::max(0., -vl0 / v); // threshold VL > 0

		const double den = v * (s + 1.) / Tp;

		if (thres >= 1)
			return 0.;

		else if (dt <= thres)
			return 0.;

		else if (dt <= 1.)
			return (pow(vl0 + v * dt, s + 1.) - pow(vl0 + v * thres, s + 1.)) / den;

		else if (dt < 1. + r * (1. - thres))
		{
			const double F1 = pow(vl0 + v, s + 1.) - pow(vl0 + v * thres, s + 1.);
			return (F1 + r * (pow(vl0 + v, s + 1) - pow(vl0 + v * (1. + r - dt) / r, s + 1))) / den;
		}
		else
		{
			const double F1 = pow(vl0 + v, s + 1.) - pow(vl0 + v * thres, s + 1.);
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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////  InfecProfileRegistry : class  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

InfecProfileRegistry::InfecProfileRegistry()
	: DistrRegistry("infectivity profile")
{
	register_distr("ferretti20", []()
				   { return std::make_shared<InfecProfile_Ferretti20>(); });
	// UPDATE: this is not a distribution...
	register_distr("cortier25", []()
				   { return std::make_shared<InfecProfile_Cortier25>(); });
}

// Singleton pattern
InfecProfileRegistry &InfecProfileRegistry::inst()
{
	static InfecProfileRegistry _inst;
	return _inst;
}
