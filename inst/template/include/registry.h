/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__REGISTRY__H
#define HOUSEPI__REGISTRY__H

#include <memory>
#include <string>
#include <vector>
#include <functional>

#include "tools/operators.h"
#include "tools/types.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  *** FunctionRegistry : class ***  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

class Function;
class Distribution;

class FunctionRegistry
{
public:
    virtual std::vector<std::string> get_options() const { return hsp::extract_keys(factory_map_); }
    virtual std::unique_ptr<Function> create(const std::string &opt, bool is_prim = false) const;

    template <typename T>
    void register_member()
    {
        factory_map_[T().name()] = []()
        { return std::make_unique<T>(); };
    }

private:
    hsp::dict<std::function<std::unique_ptr<Function>()>> factory_map_ = {};
};

class DistributionRegistry : public FunctionRegistry
{
public:
    std::vector<std::string> get_options() const override final { return hsp::extract_keys(factory_map_); }
    virtual std::unique_ptr<Distribution> create_distr(const std::string &opt, bool is_cumul = false) const;
    std::unique_ptr<Function> create(const std::string &opt, bool is_prim = false) const override final;

    template <typename T>
    void register_member()
    {
        factory_map_[T().name()] = []()
        { return std::make_unique<T>(); };
    }

private:
    hsp::dict<std::function<std::unique_ptr<Distribution>()>> factory_map_ = {};
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  *** Registry Implementations ***  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Separate functions registries.
#define DECLARE_HOUSEPI_REGISTRY(ClassName, RegistryType) \
    class ClassName : public RegistryType                 \
    {                                                     \
    public:                                               \
        static ClassName &instance()                      \
        {                                                 \
            static ClassName _inst;                       \
            return _inst;                                 \
        }                                                 \
    };
DECLARE_HOUSEPI_REGISTRY(ArithmeticRegistry, FunctionRegistry)
DECLARE_HOUSEPI_REGISTRY(TranscendentalRegistry, FunctionRegistry)
DECLARE_HOUSEPI_REGISTRY(UserFuncRegistry, FunctionRegistry)
DECLARE_HOUSEPI_REGISTRY(ContinuosDistrRegistry, DistributionRegistry)
DECLARE_HOUSEPI_REGISTRY(DiscreteDistrRegistry, DistributionRegistry)
DECLARE_HOUSEPI_REGISTRY(InfecProfileRegistry, DistributionRegistry)

/// Register a Function derived class T to FungRegistry R.
template <typename T, typename R>
struct StaticRegistatror
{
    StaticRegistatror()
    {
        R::instance().template register_member<T>();
    }
};
#define REGISTER_HOUSEPI_CLASS(T, R)       \
    namespace                              \
    {                                      \
        static StaticRegistatror<T, R>     \
            housepi_register_##T##_to_##R; \
    }

#endif // HOUSEPI__REGISTRY__H
