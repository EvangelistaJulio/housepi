/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Author: Francesco Camaglia
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 *
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 *
 * distribution_registry.h
 *
 */

#ifndef HOUSEPI__DISTRIBUTION__REGISTRY__H
#define HOUSEPI__DISTRIBUTION__REGISTRY__H

#include <memory>

#include "tools/types.h"
#include "function.h"
#include "function_registry.h"
#include "distribution.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// *** Distribution Registry : class *** ////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Distributions Registry base class.
class DistributionRegistry : public FunctionRegistry
{
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
	virtual std::unique_ptr<Distribution> get_distr(const std::string &opt, bool is_cumul = false) const;
	/// Override get method of generic function with distribution.
	std::unique_ptr<Function> get(const std::string &opt, bool is_prim = false) const override final { return get_distr(opt, is_prim); };

	/////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////  Distribution Registry : setter  ////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////

	template <typename T>
	void register_distr()
	{
		static_assert(std::is_base_of<Distribution, T>::value,
					  "T must be derived from Distribution");
		factory_map_[T().name()] = []()
		{ return std::make_unique<T>(); };
	}

private:
	std::unordered_map<std::string, std::function<std::unique_ptr<Distribution>()>> factory_map_ = {};
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
std::unique_ptr<Distribution> compose_distr_registries(const std::string &opt,
													   const std::string &name,
													   const std::vector<const DistributionRegistry *> &registries,
													   bool is_primitive = false);

#endif // HOUSEPI__DISTRIBUTION__REGISTRY__H
