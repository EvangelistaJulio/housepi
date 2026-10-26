/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "function.h"

//////////////////////////////  registry : functions and distributions  //////////////////////////////

namespace housepi
{
    class Function;
    class Distribution;
    class ModelBase;

    /// Factory-based creation of Function, Distribution and ModelBase objects by name.
    namespace registry
    {
        /// Categories a name can be looked up in.
        enum class Kind : uint8_t
        {
            Custom,
            Arithmetic,
            Transcendental,
            Continuous,
            Discrete,
            InfecProfile
        };

        /// Whether `kind` holds distributions, the others plain functions.
        constexpr bool is_distribution_kind(Kind kind)
        {
            return kind == Kind::Continuous ||
                   kind == Kind::Discrete ||
                   kind == Kind::InfecProfile;
        }

        /// Map names to factories of Function objects.
        ///
        /// * each function type registers a factory under its `NAME`;
        /// * new types are added without touching core code.
        class RegFunctions
        {
        public:
            using func_uptr = std::unique_ptr<Function>;

            /// Whether a factory is registered under `name`.
            bool contains(std::string_view name) const { return factory_map_.count(name); }

            /// Get every registered name.
            std::vector<std::string> names() const;

            /// Create the function registered under `name`, or nullptr if unknown.
            /// @param is_prim Whether the object evaluates as its primitive.
            virtual func_uptr create(std::string_view name, bool is_prim = false) const;

            /// Register the factory of function class `T` under its name.
            /// Distributions go in a distribution kind `K` and plain functions in a function kind.
            template <Kind K, typename T>
            void register_factory()
            {
                static_assert(std::is_base_of_v<Distribution, T> == is_distribution_kind(K),
                              "Class is registered under the wrong kind of registry.");
                factory_map_.emplace(std::string(T::NAME),
                                     []() -> std::unique_ptr<T>
                                     { return std::make_unique<T>(); });
            }

        protected:
            /// Factories, keyed by the registered names.
            std::map<std::string, std::function<func_uptr()>, std::less<>> factory_map_ = {};
        };

        /// Map names to factories of Distribution objects.
        class RegDistributions : public RegFunctions
        {
        public:
            /// Create the distribution registered under `name`, or nullptr if unknown.
            /// @param is_cumul Whether the object evaluates as its cumulative.
            func_uptr create(std::string_view name, bool is_cumul = false) const override final;
        };

        /// Registry holding `kind`.
        /// Defined in source, so every DLL shares one instance.
        RegFunctions &instance(Kind kind);

        /// Registries searched for a name, in order: the first carrying it wins.
        struct SearchOrder
        {
            std::vector<Kind> kinds; /// Registries, in search order.

            auto begin() const { return kinds.begin(); }
            auto end() const { return kinds.end(); }

            /// Whether some registry carries `name`.
            bool contains(std::string_view name) const;
            /// List the names registered, in order of first appearance.
            std::vector<std::string> names() const;
        };

        /// Distribution registries.
        /// * name in both registries -> the discrete one when `discr`.
        SearchOrder distr_search_order(bool discr);
        /// Deterministic function registries.
        SearchOrder func_search_order();
        /// Function, distribution and infectivity-profile registries.
        /// * name in both distribution registries -> the discrete one when `discr`.
        SearchOrder func_distr_search_order(bool discr);
        /// Infectivity-profile registry, then the distribution registries.
        /// * name in both distribution registries -> the discrete one when `discr`.
        SearchOrder infec_profile_search_order(bool discr);

        /// Create a function from the first registry carrying the name.
        /// * no registry carries it -> error suggesting the closest registered name.
        /// @param flag Evaluation flag passed to the registry: primitive for a function, cumulative for a distribution.
        std::unique_ptr<Function> create_function(
            std::string_view name,
            ErrType err_type,
            const SearchOrder &order,
            bool flag = false);

        /// Create a distribution from the first registry carrying the name.
        /// An NA name carries no definition and yields a null pointer.
        /// @param order Distribution registries only, as the distribution and profile search orders hold.
        std::unique_ptr<Distribution> create_distribution(
            std::string_view name,
            ErrType err_type,
            const SearchOrder &order,
            bool is_cumul = false);

        /// Clickable R call showing registered entry `name` and its arguments.
        /// * profile or distribution -> the browser of its registry;
        /// * any other name -> `browse_functions()`.
        std::string browse_call(std::string_view name);
    }
}

/// Register a function class into the registry of the given kind, at load time.
#define HOUSEPI_REGISTER_FUNCTION(T, kind)                                 \
    namespace                                                              \
    {                                                                      \
        struct T##kind##_AutoRegister                                      \
        {                                                                  \
            T##kind##_AutoRegister()                                       \
            {                                                              \
                housepi::registry::instance(housepi::registry::Kind::kind) \
                    .register_factory<housepi::registry::Kind::kind, T>(); \
            }                                                              \
        };                                                                 \
        static T##kind##_AutoRegister T##kind##_auto_register;             \
    }

//////////////////////////////  registry : models  //////////////////////////////

namespace housepi::registry
{
    /// Map names to factories of ModelBase objects.
    ///
    /// * model DLLs register their factory into this singleton at load time;
    /// * the singleton lives in housepi.so: every model DLL shares one instance.
    class RegModels
    {
    public:
        using model_uptr = std::unique_ptr<ModelBase>;

        /// Get the global model registry.
        /// Defined in source, so every DLL shares one instance.
        static RegModels &instance();

        /// Register the factory of model class `T` under `name`.
        template <typename T>
        void register_factory(const std::string &name)
        {
            factory_map_.emplace(name,
                                 []() -> std::unique_ptr<T>
                                 { return std::make_unique<T>(); });
        }

        /// Create the model registered under `name`.
        /// * unregistered name -> error.
        model_uptr create(std::string_view name) const;

    private:
        /// Factories, keyed by model name.
        std::map<std::string, std::function<model_uptr()>, std::less<>> factory_map_;
    };
}

/// Namespace of a module model: prefixed so any folder name yields an identifier.
#define HOUSEPI_MODEL_NAMESPACE(module) model_##module

/// Register a module model under the module name, at load time.
#define HOUSEPI_REGISTER_MODEL(module)                                                  \
    namespace                                                                           \
    {                                                                                   \
        struct model_##module##_AutoRegister                                            \
        {                                                                               \
            model_##module##_AutoRegister()                                             \
            {                                                                           \
                housepi::registry::RegModels::instance()                                \
                    .register_factory<housepi::HOUSEPI_MODEL_NAMESPACE(module)::Model>( \
                        #module);                                                       \
            }                                                                           \
        };                                                                              \
        static model_##module##_AutoRegister model_##module##_auto_register;            \
    }
