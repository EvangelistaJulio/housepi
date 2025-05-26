/*      __
 *     / /_  ____  __  __ ____ __  ____  __
 *    / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *   / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *  /_/ /_/\____/\____/____/\___/ ____/_/
 *                             /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * distribution.h
 *
 */

#ifndef HOUSEPI__DISTRIBUTION__H
#define HOUSEPI__DISTRIBUTION__H

#include <memory>
#include <string>
#include <vector>

#include "tools.h"
#include "function.h"

//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------//

#define D_DISCR_TYPE_NONE 0 // none
#define D_DISCR_TYPE_LEFT 1 // left
#define D_DISCR_TYPE_CNTR 2 // center
#define D_DISCR_TYPE_DFND 3 // defined

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  Base Distribution : class  ////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* A parent class to all distributions. */
class Distribution : public Function
{
public:
	//=======================//
	//===== CONSTRUCTOR =====//
	//=======================//

	Distribution(const std::string &name = HOUSEPI_NA_STRING,
				 const str_vec &arg_name = {},
				 const dbl_vec &arg_value = {},
				 double domain_min = HOUSEPI_NA_DOUBLE,
				 double domain_max = HOUSEPI_NA_DOUBLE,
				 bool is_min_open = true,
				 bool is_max_open = true);
	virtual ~Distribution();

	//==================//
	//===== SETTER =====//
	//==================//

	virtual void init_reshaping(bool discrete);

	//==================//
	//===== GETTER =====//
	//==================//

	bool is_discrete() const;
	double pdf(double x) const;
	double cdf(double x) const;
	double logpdf(double x) const;
	double logratio(double xold, double xnew) const;
	double rand(rnd_gen &prng);

	// note: as a function, distribution is defined by its pdf
	// FIXME: how to access cdf?
	double return_value() const override;

	//================//
	//===== SHOW =====//
	//================//

	const std::vector<dbl_vec> showcase(int n_bins = 2000,
										double ytol = 1.e-5) const;

protected:
	//==================//
	//===== SETTER =====//
	//==================//
	void init_functional_discretization(ushort discr_type);
	void init_functional_truncation(ushort discr_type);

	void init_min_zero_reshaping(bool discrete, bool is_open_zero);
	void init_min_zero_discr_defined_reshaping();

	//==================//
	//===== GETTER =====//
	//==================//
	// note: distribution specific virtual methods
	virtual double raw_pdf(double x) const;
	virtual double raw_cdf(double x) const;
	virtual double raw_logpdf(double x) const;
	virtual double raw_logratio(double xold, double xnew) const;
	virtual void update_args();
	virtual double raw_rand(rnd_gen &prng);

private:
	double _rejection_sampling(rnd_gen &prng, double (Distribution::*func)(rnd_gen &));

	//===== Functional methods for public getter methods =====//
	double (Distribution::*_do_pdf)(double) const = &Distribution::__default_pdf;
	double (Distribution::*_do_cdf)(double) const = &Distribution::__default_cdf;
	double (Distribution::*_do_logpdf)(double) const = &Distribution::__default_logpdf;
	double (Distribution::*_do_logratio)(double, double) const = &Distribution::__default_logratio;
	double (Distribution::*_do_rand)(rnd_gen &) = &Distribution::__default_rand;

	//----- Default functional methods -----//
	double __default_pdf(double x) const;
	double __default_cdf(double x) const;
	double __default_logpdf(double x) const;
	double __default_logratio(double xold, double xnew) const;
	double __default_rand(rnd_gen &prng);

	//----- Functional methods for the discrete defined distributions -----//
	double __discr_defined_pdf(double x) const;
	double __discr_defined_cdf(double x) const;
	double __discr_defined_logpdf(double x) const;
	double __discr_defined_logratio(double xold, double xnew) const;
	// note : __discr_defined_rand == __default_rand

	//----- Functional methods for the center adjusted discretization -----//
	double __discr_center_pdf(double x) const;
	double __discr_center_cdf(double x) const;
	double __discr_center_logpdf(double x) const;
	double __discr_center_logratio(double xold, double xnew) const;
	double __discr_center_rand(rnd_gen &prng);

	//----- Functional methods for the left adjusted discretization -----//
	double __discr_left_pdf(double x) const;
	double __discr_left_cdf(double x) const;
	double __discr_left_logpdf(double x) const;
	double __discr_left_logratio(double xold, double xnew) const;
	double __discr_left_rand(rnd_gen &prng);

	//===== Functional methods for truncation =====//
	double (Distribution::*_do_trunc)(double) const = &Distribution::__trunc_dummy;
	double (Distribution::*_do_logtrunc)(double) const = &Distribution::__logtrunc_dummy;
	double (Distribution::*_do_shift)(double) const = &Distribution::__shift_dummy;

	//----- Functional methods for the truncation and/or discretization -----//

	// default case
	double __trunc_dummy(double x) const;

	// continuos case
	double __trunc_min(double x) const;
	double __trunc_max(double x) const;
	double __trunc_both_sides(double x) const;
	double __trunc_rand(rnd_gen &prng);

	// case of discrete defined
	double __trunc_min_discr_defined(double x) const;
	double __trunc_max_discr_defined(double x) const;
	double __trunc_both_sides_discr_defined(double x) const;
	// note : __trunx_rand_discr_defined == __trunc_rand

	// case of center discretization
	double __trunc_min_discr_center(double x) const;
	double __trunc_max_discr_center(double x) const;
	double __trunc_both_sides_discr_center(double x) const;
	double __trunc_rand_discr_center(rnd_gen &prng);

	// case of left discretization
	double __trunc_min_discr_left(double x) const;
	double __trunc_max_discr_left(double x) const;
	double __trunc_both_sides_discr_left(double x) const;
	double __trunc_rand_discr_left(rnd_gen &prng);

	//----- Functional methods for the log-truncation -----//

	// default case
	double __logtrunc_dummy(double x) const;

	// continuos case
	double __logtrunc_min(double x) const;
	double __logtrunc_max(double x) const;
	double __logtrunc_both_sides(double x) const;

	// case of discrete defined
	double __logtrunc_min_discr_defined(double x) const;
	double __logtrunc_max_discr_defined(double x) const;
	double __logtrunc_both_sides_discr_defined(double x) const;

	// case of center discretization
	double __logtrunc_min_discr_center(double x) const;
	double __logtrunc_max_discr_center(double x) const;
	double __logtrunc_both_sides_discr_center(double x) const;

	// case of left discretization
	double __logtrunc_min_discr_left(double x) const;
	double __logtrunc_max_discr_left(double x) const;
	double __logtrunc_both_sides_discr_left(double x) const;

	//----- Functional methods to shift cumulative -----//
	double __shift_dummy(double x) const;
	double __shift_min(double x) const;
	double __shift_min_discr_defined(double x) const;
	double __shift_min_discr_center(double x) const;
	double __shift_min_discr_left(double x) const;

	//===== DEVELOPER CHECK =====//
	double __d_e_v__check__pdf__norm(int n_bins = 5000,
									 double ytol = 1.e-8) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////  DistributionRegistry : class  ///////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class DistrRegistry : public FuncRegistry
{
public:
	DistrRegistry(std::string name, str_vec options)
		: FuncRegistry(name, options) {}

	virtual std::shared_ptr<Distribution> get_distr(const std::string &opt) const = 0;
	std::shared_ptr<Function> get(const std::string &opt) const override final
	{
		return get_distr(opt);
	}
};

class ContinuosDistrRegistry : public DistrRegistry
{
public:
	ContinuosDistrRegistry()
		: DistrRegistry("continuos distribution",
						{"unif", "norm", "lnorm",
						 "beta", "gamma", "reparam_gamma",
						 "rayleigh", "weibull", "exp"}) {}

	std::shared_ptr<Distribution> get_distr(const std::string &opt) const override final;
};

class DiscreteDistrRegistry : public DistrRegistry
{
public:
	DiscreteDistrRegistry()
		: DistrRegistry("discrete distribution",
						{"pois", "binom", "nbinom", "geom"}) {}

	std::shared_ptr<Distribution> get_distr(const std::string &opt) const override final;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  InfecProfileRegistry : class  //////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class InfecProfileRegistry : public DistrRegistry
{
public:
	InfecProfileRegistry()
		: DistrRegistry("infectivity profile",
						{"ferretti20", "cortier25"}) {}

	std::shared_ptr<Distribution> get_distr(const std::string &opt) const override final;
};

class CompositeDistrRegistry : public DistrRegistry
{
private:
	std::vector<std::shared_ptr<DistrRegistry>> _registries;

public:
	CompositeDistrRegistry(std::string name,
						   std::vector<std::shared_ptr<DistrRegistry>> registries)
		: DistrRegistry(name, {}), _registries(std::move(registries))
	{
		for (const auto &r : _registries)
		{
			const str_vec &this_opts = r->get_options();
			_options.insert(_options.end(), this_opts.begin(), this_opts.end());
		}
	}

	std::shared_ptr<Distribution> get_distr(const std::string &opt) const override final
	{
		if (tools::is_na(opt))
			return std::make_shared<Distribution>();
		else
			for (const auto &r : _registries)
				if (r->valid_option(opt))
					return r->get_distr(opt);

		throw std::invalid_argument(
			"Unknown " + _name + " \"" + opt + "\"." + tools::find_closest(opt, _options));
		return nullptr;
	}
};

#endif // HOUSEPI__DISTRIBUTION__H
