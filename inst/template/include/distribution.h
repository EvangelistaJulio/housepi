/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
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

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////// *** Distribution : class *** ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Parent class of all distributions.
class Distribution : public Function
{
public:
	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  Distribution : constructor  ////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	Distribution(const std::string &name = HSP_NA_STRING,
				 const str_vec &arg_name = {}, const dbl_vec &arg_value = {},
				 double domain_min = HSP_NA_DOUBLE, bool is_min_open = true,
				 double domain_max = HSP_NA_DOUBLE, bool is_max_open = true)
		: Function(name, arg_name, arg_value,
				   domain_min, is_min_open,
				   domain_max, is_max_open)
	{
		err_type_ = "distribution";
		set_eval_as_distr_call();
	}
	virtual ~Distribution() = default;

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////  Distribution : getter  //////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double rand(rnd_gen &prng);
	/// Return value of probability density/mass function.
	double pdf(double x) const { return eval_as_distr_call_(x); }
	/// Return value of log-probability density/mass function.
	double logpdf(double x) const { return (this->*log_func_)(x); }
	/// Return value of cumulative distribution.
	double cdf(double x) const { return eval_as_distr_primitive_(x); }
	/// Return value of log-ratio of probability density/mass functions.
	double logratio(double xo, double xn) const { return (xo != xn) ? (this->*log_ratio_)(xo, xn) : 0.; }

protected:
	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Distribution : setter [protected] /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Initialize functional discretization.
	void init_discr_helper_(unsigned short discr_type) override final;
	void init_trunc_helper_();

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Distribution : virtual [protected] ////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	virtual void update_args_() {}; // UPDATE

	/// Raw random sampling method of the distribution.
	/// @param prng : pseudo-random number generator.
	virtual double raw_rand_(rnd_gen &) { throw std::invalid_argument(
		"Wrong call to undefined method \"rand\" of " + err_name() + "."); }

	/// Raw probability density/mass function of the distribution.
	virtual double raw_pdf_(double) const { throw std::invalid_argument(
		"Wrong call to undefined method \"pdf\" of " + err_name() + "."); }

	/// Raw log-probability density/mass function of the distribution.
	virtual double raw_logpdf_(double) const { throw std::invalid_argument(
		"Wrong call to undefined method \"logpdf\" of " + err_name() + "."); }

	/// Raw cumulative distribution function.
	virtual double raw_cdf_(double) const { throw std::invalid_argument(
		"Wrong call to undefined method \"cdf\" of " + err_name() + "."); }

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////  Distribution : getter [protected] /////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Rejection sampling.
	/// @param prng : pseudo-random number generator.
	/// @param func : pointer to the raw random sampling method of the distribution.
	double rejec_sample_(rnd_gen &prng, double (Distribution::*func)(rnd_gen &));

	double (Distribution::*rand_)(rnd_gen &) = &Distribution::raw_rand_;
	double do_discr_center_rand_(rnd_gen &prng) { return static_cast<double>(std::lround(raw_rand_(prng))); }
	double do_discr_left_rand_(rnd_gen &prng) { return std::floor(raw_rand_(prng)); }

	double raw_function_(double x) const override final { return (this->*trunc_)(raw_pdf_(x)); }
	double do_function_discr_center_(double x) const override final
	{
		const double i = std::round(x);
		return (this->*trunc_)(raw_cdf_(i + .5) - raw_cdf_(i - .5));
	}
	double do_function_discr_left_(double x) const override final
	{
		const double i = std::floor(x);
		return (this->*trunc_)(raw_cdf_(i + 1.) - raw_cdf_(i));
	}

	double raw_logfunc_(double x) const override final { return (this->*logtrunc_)(raw_logpdf_(x)); }

	double raw_primitive_(double x) const override final { return (this->*trunc_)((this->*shift_)(raw_cdf_(x))); }
	double do_primitive_discr_center_(double x) const override final { return (this->*trunc_)((this->*shift_)(raw_cdf_(std::round(x) + .5))); }
	double do_primitive_discr_left_(double x) const override final { return (this->*trunc_)((this->*shift_)(raw_cdf_(std::floor(x) + 1.))); }

	double do_logratio_discr_center_(double xo, double xn) const override final
	{
		const double io = std::round(xo), in = std::round(xn);
		return std::log((raw_cdf_(in + .5) - raw_cdf_(in - .5)) / (raw_cdf_(io + .5) - raw_cdf_(io - .5)));
	}
	double do_logratio_discr_left_(double xo, double xn) const override final
	{
		const double io = std::floor(xo), in = std::floor(xn);
		return std::log((raw_cdf_(in + 1.) - raw_cdf_(in)) / (raw_cdf_(io + 1.) - raw_cdf_(io)));
	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////  Distribution : functional truncation [private] //////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double (Distribution::*trunc_)(double) const = &Distribution::do_trunc_dummy_;

	// default case
	double do_trunc_dummy_(double x) const { return x; }

	// continuos case
	double do_trunc_min_(double x) const { return x / (1. - raw_cdf_(get_min())); }
	double do_trunc_max_(double x) const { return x / (raw_cdf_(get_max())); }
	double do_trunc_both_sides_(double x) const { return x / (raw_cdf_(get_max()) - raw_cdf_(get_min())); }
	double do_trunc_rand_(rnd_gen &prng) { return rejec_sample_(prng, &Distribution::raw_rand_); }

	// case of discrete defined
	double do_trunc_min_discr_defined_(double x) const { return x / (1. - raw_cdf_(get_min() - 1)); }
	double do_trunc_max_discr_defined_(double x) const { return x / (raw_cdf_(get_max() + 1)); }
	double do_trunc_both_sides_discr_defined_(double x) const { return x / (raw_cdf_(get_max() + 1) - raw_cdf_(get_min() - 1)); }
	// note : do_trunx_rand_discr_defined == do_trunc_rand_

	// case of center discretization
	double do_trunc_min_discr_center_(double x) const { return x / (1. - raw_cdf_(std::ceil(get_min()) - .5)); }
	double do_trunc_max_discr_center_(double x) const { return x / (raw_cdf_(std::floor(get_max()) + .5)); }
	double do_trunc_both_sides_discr_center_(double x) const { return x / (raw_cdf_(std::floor(get_max()) + .5) - raw_cdf_(std::ceil(get_min()) - .5)); }
	double do_trunc_rand_discr_center_(rnd_gen &prng) { return rejec_sample_(prng, &Distribution::do_discr_center_rand_); }

	// case of left discretization
	double do_trunc_min_discr_left_(double x) const { return x / (1. - raw_cdf_(std::floor(get_min()))); }
	double do_trunc_max_discr_left_(double x) const { return x / (raw_cdf_(std::ceil(get_max() + 1))); }
	double do_trunc_both_sides_discr_left_(double x) const { return x / (raw_cdf_(std::ceil(get_max() + 1)) - raw_cdf_(std::floor(get_min()))); }
	double do_trunc_rand_discr_left_(rnd_gen &prng) { return rejec_sample_(prng, &Distribution::do_discr_left_rand_); }

	/////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////  Distribution : functional log-truncation [private] ////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// note : this is not improving the case of discretization with repspect to simply log of pdf
	double (Distribution::*logtrunc_)(double) const = &Distribution::do_logtrunc_dummy_;

	// default case
	double do_logtrunc_dummy_(double x) const { return x; }

	// continuos case
	double do_logtrunc_min_(double x) const { return x - std::log(1. - raw_cdf_(get_min())); }
	double do_logtrunc_max_(double x) const { return x - std::log(raw_cdf_(get_max())); }
	double do_logtrunc_both_sides_(double x) const { return x - std::log(raw_cdf_(get_max()) - raw_cdf_(get_min())); }

	// case of discrete defined
	double do_logtrunc_min_discr_defined_(double x) const { return x - std::log(1. - raw_cdf_(get_min() - 1)); }
	double do_logtrunc_max_discr_defined_(double x) const { return x - std::log(raw_cdf_(get_max() + 1)); }
	double do_logtrunc_both_sides_discr_defined_(double x) const { return x - std::log(raw_cdf_(get_max() + 1) - raw_cdf_(get_min() - 1)); }

	// case of center discretization
	double do_logtrunc_min_discr_center_(double x) const { return x - std::log(1. - raw_cdf_(std::ceil(get_min()) - .5)); }
	double do_logtrunc_max_discr_center_(double x) const { return x - std::log(raw_cdf_(std::floor(get_max()) + .5)); }
	double do_logtrunc_both_sides_discr_center_(double x) const { return x - std::log(raw_cdf_(std::floor(get_max()) + .5) - raw_cdf_(std::ceil(get_min()) - .5)); }

	// case of left discretization
	double do_logtrunc_min_discr_left_(double x) const { return x - std::log(1. - raw_cdf_(std::floor(get_min()))); }
	double do_logtrunc_max_discr_left_(double x) const { return x - std::log(raw_cdf_(std::ceil(get_max()))); }
	double do_logtrunc_both_sides_discr_left_(double x) const { return x - std::log(raw_cdf_(std::ceil(get_max())) - raw_cdf_(std::floor(get_min()))); }

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////  Distribution : functional shift [private] ////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	double (Distribution::*shift_)(double) const = &Distribution::do_shift_dummy_;
	double do_shift_dummy_(double x) const { return x; }
	double do_shift_min_(double x) const { return x - raw_cdf_(get_min()); }
	double do_shift_min_discr_defined_(double x) const { return x - raw_cdf_(get_min() - 1); }
	double do_shift_min_discr_center_(double x) const { return x - raw_cdf_(std::floor(get_min()) - .5); }
	double do_shift_min_discr_left_(double x) const { return x - raw_cdf_(std::floor(get_min())); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// *** Distribution Registry : class *** ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Distributions Registry base class.
class DistributionRegistry : public FunctionRegistry
{
	using FactoryFunc = std::function<std::shared_ptr<Function>()>;
	using FactoryDistr = std::function<std::shared_ptr<Distribution>()>;

public:
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////  Distribution Registry : constructor  /////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	DistributionRegistry(const std::string &name)
		: FunctionRegistry(name) {}
	virtual ~DistributionRegistry() = default;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  Distribution Registry : getter  ////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	/// Override get_options method of generic function with distribution.
	str_vec get_options() const override final;
	/// Method to get a distribution.
	virtual std::shared_ptr<Distribution> get_distr(const std::string &opt) const;
	/// Override get method of generic function with distribution.
	std::shared_ptr<Function> get(const std::string &opt, bool is_primitive) const override final;

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  Distribution Registry : setter  ////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	void register_distr()
	{
		static_assert(std::is_base_of<Distribution, T>::value,
					  "T must be derived from Distribution");
		factory_map_[T().name()] = []()
		{ return std::make_shared<T>(); };
	}

private:
	std::unordered_map<std::string, FactoryDistr> factory_map_ = {};
};

/// Register a Distribution derived class T to DistributionRegistry derived class C.
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

/// Register a Distribution derived class T to DistributionRegistry derived class R.
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

/// Compose registries to get a distribution from them.
std::shared_ptr<Distribution> compose_distr_registries(const std::string &opt,
													   const std::string &name,
													   const std::vector<const DistributionRegistry *> &registries);

#endif // HOUSEPI__DISTRIBUTION__H
