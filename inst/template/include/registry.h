/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

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
    using func_uptr = std::unique_ptr<Function>;

    size_t size() const noexcept { return factory_map_.size(); }
    bool count(std::string_view opt) const { return factory_map_.count(opt); }
    std::vector<std::string> get_options() const;
    virtual func_uptr create(std::string_view opt, bool is_prim = false) const;

    /// UPDATE: this requires to create an instance of T.
    template <typename T>
    void register_member()
    {
        factory_map_[std::string(T().name())] = []()
        { return std::make_unique<T>(); };
    }

protected:
    std::map<std::string, std::function<func_uptr()>, std::less<>> factory_map_ = {};
};

class DistributionRegistry : public FunctionRegistry
{
public:
    func_uptr create(std::string_view opt, bool is_prim = false) const override final;
};

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  *** Registry Implementations ***  /////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

/// Declare a housepi functions registry.
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

#define REGISTER_HOUSEPI_FUNCTION(T, R)                                \
    namespace                                                          \
    {                                                                  \
        struct T##_##R##_AutoRegister                                  \
        {                                                              \
            T##_##R##_AutoRegister()                                   \
            {                                                          \
                R::instance().template register_member<T>();           \
            }                                                          \
        };                                                             \
        static T##_##R##_AutoRegister T##_##R##_housepi_auto_register; \
    }
