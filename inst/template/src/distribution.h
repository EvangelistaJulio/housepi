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

	Distribution(const std::string &name = HSP_NA_STRING,
				 const str_vec &arg_name = {},
				 const dbl_vec &arg_value = {},
				 double domain_min = HSP_NA_DOUBLE,
				 double domain_max = HSP_NA_DOUBLE,
				 bool is_min_open = true,
				 bool is_max_open = true)
		: Function(name, arg_name, arg_value,
				   domain_min, domain_max, is_min_open, is_max_open) {}
	virtual ~Distribution() = default;

	//==================//
	//===== GETTER =====//
	//==================//

	bool is_discrete() const { return pdf_ != &Distribution::do_default_pdf_; }
	double logpdf(double x) const { return (this->*logpdf_)(x); }
	double logratio(double xo, double xn) const { return (xo != xn) ? (this->*logratio_)(xo, xn) : 0.; }
	double rand(rnd_gen &prng);

	/// Return value of probability density
	double pdf(double x) const { return is_valid(x) ? (this->*pdf_)(x) : 0.; }

	/// Return value of cumulative distribution
	double cdf(double x) const { return is_valid(x) ? (this->*cdf_)(x) : (is_above(x) ? (this->*cdf_)(get_max()) : 0.); }

	// note: as a function, distribution is defined by its pdf or cdf
	double return_value(double x) const override { return (this->*return_value_)(x); }

	//==================//
	//===== SETTER =====//
	//==================//

	void set_return_cumul() { return_value_ = &Distribution::cdf; }
	virtual void init_reshaping(bool discrete) override;

protected:
	//==================//
	//===== SETTER =====//
	//==================//
	void init_functional_discretization(unsigned short discr_type);
	void init_functional_truncation(unsigned short discr_type);

	bool init_range_min(double min_val, bool is_open_zero);
	void init_reshaping_min_zero(bool discrete, bool is_open_zero);
	void init_reshaping_min_zero_discr_defined();

	//==================//
	//===== GETTER =====//
	//==================//

	virtual void update_args() {};

	virtual double raw_rand(rnd_gen &) { throw std::invalid_argument(
		"Called undefined method \"rand\" of distribution \"" + name() + "\"."); }

	virtual double raw_pdf(double) const { throw std::invalid_argument(
		"Called undefined method \"pdf\" of distribution \"" + name() + "\"."); }

	virtual double raw_cdf(double) const { throw std::invalid_argument(
		"Called undefined method \"cdf\" of distribution \"" + name() + "\"."); }

	virtual double raw_logpdf(double) const { throw std::invalid_argument(
		"Called undefined method \"logpdf\" of distribution \"" + name() + "\"."); }

	virtual double raw_logratio(double, double) const { return 0.; }

private:
	double rejec_sample_(rnd_gen &prng, double (Distribution::*func)(rnd_gen &));

	//===== Functional methods for return_value =====//
	double (Distribution::*return_value_)(double) const = &Distribution::pdf;

	//===== Functional methods for public getter methods =====//
	double (Distribution::*pdf_)(double) const = &Distribution::do_default_pdf_;
	double (Distribution::*cdf_)(double) const = &Distribution::do_default_cdf_;
	double (Distribution::*logpdf_)(double) const = &Distribution::do_default_logpdf_;
	double (Distribution::*logratio_)(double, double) const = &Distribution::do_default_logratio_;
	double (Distribution::*rand_)(rnd_gen &) = &Distribution::do_default_rand_;

	//===== Default functional methods =====//
	double do_default_pdf_(double x) const { return (this->*trunc_)(raw_pdf(x)); }
	double do_default_cdf_(double x) const { return (this->*trunc_)((this->*shift_)(raw_cdf(x))); }
	double do_default_logpdf_(double x) const { return (this->*logtrunc_)(raw_logpdf(x)); }
	double do_default_logratio_(double xo, double xn) const { return raw_logratio(xo, xn); }
	double do_default_rand_(rnd_gen &prng) { return raw_rand(prng); }

	//===== Functional methods for the discrete defined distributions =====//
	double do_discr_defined_pdf_(double x) const { return do_default_pdf_(std::floor(x)); }
	double do_discr_defined_cdf_(double x) const { return do_default_cdf_(std::floor(x)); }
	double do_discr_defined_logpdf_(double x) const { return do_default_logpdf_(std::floor(x)); }
	double do_discr_defined_logratio_(double xo, double xn) const { return do_default_logratio_(std::floor(xo), std::floor(xn)); }
	// note : do_discr_defined_rand == do_default_rand_

	//===== Functional methods for the center adjusted discretization =====//
	double do_discr_center_pdf_(double x) const
	{
		const double i = std::round(x);
		return (this->*trunc_)(raw_cdf(i + .5) - raw_cdf(i - .5));
	}
	double do_discr_center_cdf_(double x) const { return (this->*trunc_)((this->*shift_)(raw_cdf(std::round(x) + .5))); }
	double do_discr_center_logpdf_(double x) const { return std::log(do_discr_center_pdf_(x)); }
	double do_discr_center_logratio_(double xo, double xn) const
	{
		const double io = std::round(xo), in = std::round(xn);
		return std::log((raw_cdf(in + .5) - raw_cdf(in - .5)) / (raw_cdf(io + .5) - raw_cdf(io - .5)));
	}
	double do_discr_center_rand_(rnd_gen &prng) { return static_cast<double>(std::lround(do_default_rand_(prng))); }

	//===== Functional methods for the left adjusted discretization =====//
	double do_discr_left_pdf_(double x) const
	{
		const double i = std::floor(x);
		return (this->*trunc_)(raw_cdf(i + 1.) - raw_cdf(i));
	}
	double do_discr_left_cdf_(double x) const { return (this->*trunc_)((this->*shift_)(raw_cdf(std::floor(x) + 1.))); }
	double do_discr_left_logpdf_(double x) const { return std::log(do_discr_left_pdf_(x)); }
	double do_discr_left_logratio_(double xo, double xn) const
	{
		const double io = std::floor(xo), in = std::floor(xn);
		return std::log((raw_cdf(in + 1.) - raw_cdf(in)) / (raw_cdf(io + 1.) - raw_cdf(io)));
	}
	double do_discr_left_rand_(rnd_gen &prng) { return std::floor(do_default_rand_(prng)); }

	//===== Functional methods for truncation =====//
	double (Distribution::*trunc_)(double) const = &Distribution::do_trunc_dummy_;
	double (Distribution::*logtrunc_)(double) const = &Distribution::do_logtrunc_dummy_;
	double (Distribution::*shift_)(double) const = &Distribution::do_shift_dummy_;

	//----- Functional methods for the truncation and/or discretization -----//

	// default case
	double do_trunc_dummy_(double x) const { return x; }

	// continuos case
	double do_trunc_min_(double x) const { return x / (1. - raw_cdf(get_min())); }
	double do_trunc_max_(double x) const { return x / (raw_cdf(get_max())); }
	double do_trunc_both_sides_(double x) const { return x / (raw_cdf(get_max()) - raw_cdf(get_min())); }
	double do_trunc_rand_(rnd_gen &prng) { return rejec_sample_(prng, &Distribution::do_default_rand_); }

	// case of discrete defined
	double do__trunc_min_discr_defined_(double x) const { return x / (1. - raw_cdf(get_min() - 1)); }
	double do_trunc_max_discr_defined_(double x) const { return x / (raw_cdf(get_max() + 1)); }
	double do_trunc_both_sides_discr_defined_(double x) const { return x / (raw_cdf(get_max() + 1) - raw_cdf(get_min() - 1)); }
	// note : do_trunx_rand_discr_defined == do_trunc_rand_

	// case of center discretization
	double do_trunc_min_discr_center_(double x) const { return x / (1. - raw_cdf(std::ceil(get_min()) - .5)); }
	double do_trunc_max_discr_center_(double x) const { return x / (raw_cdf(std::floor(get_max()) + .5)); }
	double do_trunc_both_sides_discr_center_(double x) const { return x / (raw_cdf(std::floor(get_max()) + .5) - raw_cdf(std::ceil(get_min()) - .5)); }
	double do_trunc_rand_discr_center_(rnd_gen &prng) { return rejec_sample_(prng, &Distribution::do_discr_center_rand_); }

	// case of left discretization
	double do_trunc_min_discr_left_(double x) const { return x / (1. - raw_cdf(std::floor(get_min()))); }
	double do_trunc_max_discr_left_(double x) const { return x / (raw_cdf(std::ceil(get_max() + 1))); }
	double do_trunc_both_sides_discr_left_(double x) const { return x / (raw_cdf(std::ceil(get_max() + 1)) - raw_cdf(std::floor(get_min()))); }
	double do_trunc_rand_discr_left_(rnd_gen &prng) { return rejec_sample_(prng, &Distribution::do_discr_left_rand_); }

	//----- Functional methods for the log-truncation -----//

	// default case
	double do_logtrunc_dummy_(double x) const { return x; }

	// continuos case
	double do_logtrunc_min_(double x) const { return x - std::log(1. - raw_cdf(get_min())); }
	double do_logtrunc_max_(double x) const { return x - std::log(raw_cdf(get_max())); }
	double do_logtrunc_both_sides_(double x) const { return x - std::log(raw_cdf(get_max()) - raw_cdf(get_min())); }

	// case of discrete defined
	double do_logtrunc_min_discr_defined_(double x) const { return x - std::log(1. - raw_cdf(get_min() - 1)); }
	double do_logtrunc_max_discr_defined_(double x) const { return x - std::log(raw_cdf(get_max() + 1)); }
	double do_logtrunc_both_sides_discr_defined_(double x) const { return x - std::log(raw_cdf(get_max() + 1) - raw_cdf(get_min() - 1)); }

	// case of center discretization
	double do_logtrunc_min_discr_center_(double x) const { return x - std::log(1. - raw_cdf(std::ceil(get_min()) - .5)); }
	double do_logtrunc_max_discr_center_(double x) const { return x - std::log(raw_cdf(std::floor(get_max()) + .5)); }
	double do_logtrunc_both_sides_discr_center_(double x) const { return x - std::log(raw_cdf(std::floor(get_max()) + .5) - raw_cdf(std::ceil(get_min()) - .5)); }

	// case of left discretization
	double do_logtrunc_min_discr_left_(double x) const { return x - std::log(1. - raw_cdf(std::floor(get_min()))); }
	double do_logtrunc_max_discr_left_(double x) const { return x - std::log(raw_cdf(std::ceil(get_max()))); }
	double do_logtrunc_both_sides_discr_left_(double x) const { return x - std::log(raw_cdf(std::ceil(get_max())) - raw_cdf(std::floor(get_min()))); }

	//----- Functional methods to shift cumulative -----//
	double do_shift_dummy_(double x) const { return x; }
	double do_shift_min_(double x) const { return x - raw_cdf(get_min()); }
	double do_shift_min_discr_defined_(double x) const { return x - raw_cdf(get_min() - 1); }
	double do_shift_min_discr_center_(double x) const { return x - raw_cdf(std::floor(get_min()) - .5); }
	double do_shift_min_discr_left_(double x) const { return x - raw_cdf(std::floor(get_min())); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////  Distribution Registries  //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/* Distributions Registry base class. */
class DistributionRegistry : public FunctionRegistry
{
	using FactoryFunc = std::function<std::shared_ptr<Function>()>;
	using FactoryDistr = std::function<std::shared_ptr<Distribution>()>;

public:
	DistributionRegistry(const std::string &name)
		: FunctionRegistry(name) {}
	virtual ~DistributionRegistry() = default;

	//==================//
	//===== GETTER =====//
	//==================//

	/* Override get_options method of generic function with distribution. */
	str_vec get_options() const override final;

	/* Override get method of generic function with distribution. */
	std::shared_ptr<Function> get(const std::string &opt) const override final;

	/* Override valid options to allow for cumulative distribution call.*/
	bool is_valid_option(const std::string &opt) const override final;

	virtual std::shared_ptr<Distribution> get_distr(const std::string &opt) const;

	//==================//
	//===== SETTER =====//
	//==================//

	template <typename T>
	void register_distr()
	{
		static_assert(std::is_base_of<Distribution, T>::value,
					  "T must be derived from Distribution");
		_factory_map[T().name()] = []()
		{ return std::make_shared<T>(); };
	}

private:
	std::unordered_map<std::string, FactoryDistr> _factory_map = {};
};

/* Register a Distribution derived class T to DistributionRegistry derived class C. */
template <typename C, typename T>
void register_distr()
{
	static_assert(std::is_base_of<DistributionRegistry, C>::value,
				  "C must be derived from DistributionRegistry");
	C::inst().template register_distr<T>();
}

#define DECLARE_HOUSEPI_DISTRIBUTION_REGISTRY(ClassName, Description) \
	class ClassName : public DistributionRegistry                     \
	{                                                                 \
	public:                                                           \
		static ClassName &inst();                                     \
                                                                      \
	private:                                                          \
		ClassName() : DistributionRegistry(Description) {}            \
	};
DECLARE_HOUSEPI_DISTRIBUTION_REGISTRY(ContinuosDistrRegistry, "continuos distribution")
DECLARE_HOUSEPI_DISTRIBUTION_REGISTRY(DiscreteDistrRegistry, "discrete distribution")
DECLARE_HOUSEPI_DISTRIBUTION_REGISTRY(InfecProfileRegistry, "infectivity profile")

/* Register a Distribution derived class T to DistributionRegistry derived class R. */
template <typename R, typename T>
struct StaticDistributionRegistatror
{
	StaticDistributionRegistatror()
	{
		static_assert(std::is_base_of<DistributionRegistry, R>::value,
					  "R must be derived from DistributionRegistry.");
		R::inst().template register_distr<T>();
	}
};
#define REGISTER_HOUSEPI_DISTRIBUTION(R, T)        \
	namespace                                      \
	{                                              \
		static StaticDistributionRegistatror<R, T> \
			housepi_register_##T##_to_##R;         \
	}
#define REGISTER_USER_CONTINUOS_DISTRIBUTION(T) REGISTER_HOUSEPI_DISTRIBUTION(ContinuosDistrRegistry, T)
#define REGISTER_USER_DISCRETE_DISTRIBUTION(T) REGISTER_HOUSEPI_DISTRIBUTION(DiscreteDistrRegistry, T)
#define REGISTER_USER_INFEC_PROFILE(T) REGISTER_HOUSEPI_DISTRIBUTION(InfecProfileRegistry, T)

/* Compose registries to get a distribution from them. */
std::shared_ptr<Distribution> compose_distr_registries(const std::string &opt,
													   const std::string &name,
													   const std::vector<const DistributionRegistry *> &registries);

#endif // HOUSEPI__DISTRIBUTION__H
