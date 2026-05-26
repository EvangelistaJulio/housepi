/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
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

namespace housepi
{
    class Function;
    class Distribution;

    /// @class FunctionRegistry
    /// @brief Base registry class for managing factory-based creation of Function objects.
    ///
    /// FunctionRegistry implements the registry pattern for polymorphic object creation.
    /// Derived classes maintain a factory map that enables dynamic instantiation of registered
    /// function types by name. This pattern supports extensibility without modifying core code.
    ///
    /// @details
    /// - Uses a std::map with std::function wrappers to store factory lambdas
    /// - Provides template-based registration mechanism for derived classes
    /// - Supports lookup and creation of objects by string identifier
    class FunctionRegistry
    {
    public:
        using func_uptr = std::unique_ptr<Function>;

        /// @brief Check if a factory for the given option exists.
        /// @param opt The option name to check.
        /// @return true if the option is registered, false otherwise.
        bool count(std::string_view opt) const { return factory_map_.count(opt); }

        /// @brief Get all registered option names.
        /// @return A vector of registered option identifiers.
        std::vector<std::string> get_options() const;

        /// @brief Create a Function instance for the given option.
        /// @param opt The option/type name to instantiate.
        /// @param is_prim Whether the object should be created in primitive mode (default: false).
        /// @return A unique_ptr to the created Function object.
        /// @throws std::exception if the option is not registered.
        virtual func_uptr create(std::string_view opt, bool is_prim = false) const;

        /// @brief Register a function type with the factory.
        /// @tparam T The function class type to register. Must have a name() method.
        /// @details This template method automatically creates and stores a factory lambda
        /// that can be invoked later to create instances of type T.
        template <typename T>
        void register_factory()
        {
            factory_map_.emplace(std::string(T().name()),
                                 []() -> std::unique_ptr<T>
                                 { return std::make_unique<T>(); });
        }

    protected:
        /// @brief Factory map storing function creation lambdas.
        /// @details Keys are the string names of registered types, values are factory functions
        /// that return unique_ptr<Function> instances.
        std::map<std::string, std::function<func_uptr()>, std::less<>> factory_map_ = {};
    };

    /// @class DistributionRegistry
    /// @brief Specialized registry for creating Distribution objects.
    ///
    /// DistributionRegistry extends FunctionRegistry to provide specialized factory creation
    /// for Distribution types. It overrides the create() method with final to enforce specific
    /// distribution creation semantics in derived registries.
    class DistributionRegistry : public FunctionRegistry
    {
    public:
        /// @brief Create a Distribution instance for the given option.
        /// @param opt The distribution type name to instantiate.
        /// @param is_prim Whether the object should be created in primitive mode (default: false).
        /// @return A unique_ptr to the created Function object (polymorphic base).
        /// @details This override is marked final to prevent further specialization in derived
        /// distribution registry classes.
        func_uptr create(std::string_view opt, bool is_prim = false) const override final;
    };
}

#define HOUSEPI_OVERRIDE_NAME_METHOD(name_) \
    std::string_view name() const noexcept override final { return #name_; }
#define HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_) \
    std::unique_ptr<housepi::Function> clone() const override final { return std::make_unique<Class_>(*this); }
#define HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Class_) \
    std::unique_ptr<housepi::Distribution> clone_distrib() const override final { return std::make_unique<Class_>(*this); }

/// Declare housepi functions registry.
#define HOUSEPI_DECLARE_REGISTRY_(ClassName, Type)             \
    class ClassName##Registry : public housepi::Type##Registry \
    {                                                          \
    public:                                                    \
        static ClassName##Registry &instance();                \
    };

/// Define instance method of registry.
/// NOTE: I define instance in source to avoid duplication across DLLs.
#define HOUSEPI_DEFINE_REGISTRY_INSTANCE_(ClassName)     \
    ClassName##Registry &ClassName##Registry::instance() \
    {                                                    \
        static ClassName##Registry _inst;                \
        return _inst;                                    \
    }

HOUSEPI_DECLARE_REGISTRY_(Custom, Function)
HOUSEPI_DECLARE_REGISTRY_(Arithmetic, Function)
HOUSEPI_DECLARE_REGISTRY_(Transcendental, Function)
HOUSEPI_DECLARE_REGISTRY_(ContinuosDistrib, Distribution)
HOUSEPI_DECLARE_REGISTRY_(DiscreteDistrib, Distribution)
HOUSEPI_DECLARE_REGISTRY_(InfecProfile, Distribution)

#define HOUSEPI_REGISTER_FUNCTION(T, R)                                 \
    namespace                                                           \
    {                                                                   \
        struct T##R##_AutoRegister                                      \
        {                                                               \
            T##R##_AutoRegister()                                       \
            {                                                           \
                R##Registry::instance().template register_factory<T>(); \
            }                                                           \
        };                                                              \
        static T##R##_AutoRegister T##R##_auto_register;                \
    }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////  *** ModelRegistry : class ***  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    class ModelBase;

    /// @class ModelRegistry
    /// @brief Registry for dynamic model creation by name.
    ///
    /// Model DLLs auto-register their factory into this singleton at load time.
    /// The singleton lives in housepi.so, so all model DLLs share the same instance.
    class ModelRegistry
    {
    public:
        using func_uptr = std::unique_ptr<ModelBase>;

        /// @brief Get the global ModelRegistry instance.
        /// NOTE: I define instance in source to avoid duplication across DLLs.
        static ModelRegistry &instance();

        /// @brief Register a model factory under the given name.
        /// @tparam T The model class type to register.
        /// @param name The name of the model.
        template <typename T>
        void register_factory(const std::string &name)
        {
            factory_map_.emplace(name,
                                 []() -> std::unique_ptr<T>
                                 { return std::make_unique<T>(); });
        }

        /// @brief Create a model by name.
        /// @throws std::runtime_error if the name is not registered.
        func_uptr create(std::string_view name) const;

    private:
        std::map<std::string, std::function<func_uptr()>, std::less<>> factory_map_;
    };
}

#define HOUSEPI_REGISTER_MODEL(module)                                                        \
    namespace                                                                                 \
    {                                                                                         \
        struct module##_ModelAutoRegister                                                     \
        {                                                                                     \
            module##_ModelAutoRegister()                                                      \
            {                                                                                 \
                housepi::ModelRegistry::instance().register_factory<housepi::module::Model>(#module); \
            }                                                                                 \
        };                                                                                    \
        static module##_ModelAutoRegister module##_model_auto_register;                       \
    }
