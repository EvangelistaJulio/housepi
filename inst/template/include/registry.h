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

    virtual std::string_view name() const = 0;
    bool count(std::string_view opt) const { return factory_map_.count(opt); }
    std::vector<std::string> get_options() const;
    virtual func_uptr create(std::string_view opt, bool is_prim = false) const;

    /// UPDATE: this requires to create an instance of T.
    template <typename T>
    void register_member()
    {
        factory_map_.emplace(std::string(T().name()), []()
                             { return std::make_unique<T>(); });
    }

protected:
    std::map<std::string, std::function<func_uptr()>, std::less<>> factory_map_ = {};
};

class DistributionRegistry : public FunctionRegistry
{
public:
    func_uptr create(std::string_view opt, bool is_prim = false) const override final;
};

#define HOUSEPI_OVERRIDE_NAME_METHOD(name_) \
    std::string_view name() const noexcept override final { return #name_; }
#define HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_) \
    std::unique_ptr<Function> clone() const override final { return std::make_unique<Class_>(*this); }
#define HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Class_) \
    std::unique_ptr<Distribution> clone_distrib() const override final { return std::make_unique<Class_>(*this); }

/// Declare housepi functions registry.
#define HOUSEPI_DECLARE_REGISTRY_(ClassName, Type)    \
    class ClassName##Registry : public Type##Registry \
    {                                                 \
    public:                                           \
        HOUSEPI_OVERRIDE_NAME_METHOD(ClassName)       \
        static ClassName##Registry &instance();       \
    };

/// Define instance method of registry.
#define HOUSEPI_DEFINE_REGISTRY_INSTANCE_(ClassName)     \
    ClassName##Registry &ClassName##Registry::instance() \
    {                                                    \
        static ClassName##Registry _inst;                \
        return _inst;                                    \
    }

#define HOUSEPI_REGISTER_FUNCTION(T, R)                                \
    namespace                                                          \
    {                                                                  \
        struct T##R##_AutoRegister                                     \
        {                                                              \
            T##R##_AutoRegister()                                      \
            {                                                          \
                R##Registry::instance().template register_member<T>(); \
            }                                                          \
        };                                                             \
        static T##R##_AutoRegister T##R##_auto_register;               \
    }

HOUSEPI_DECLARE_REGISTRY_(Custom, Function)
HOUSEPI_DECLARE_REGISTRY_(Arithmetic, Function)
HOUSEPI_DECLARE_REGISTRY_(Transcendental, Function)
HOUSEPI_DECLARE_REGISTRY_(ContinuosDistrib, Distribution)
HOUSEPI_DECLARE_REGISTRY_(DiscreteDistrib, Distribution)
HOUSEPI_DECLARE_REGISTRY_(InfecProfile, Distribution)
