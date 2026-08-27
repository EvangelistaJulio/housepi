/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "function.h"

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  *** registry : functions and distributions ***  ///////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi
{
    class Function;
    class Distribution;
    class ModelBase;

    /// Factory-based creation of Function, Distribution and ModelBase objects by name.
    namespace registry
    {
        /// @class RegFunctions
        /// @brief Name to factory map for Function objects.
        ///
        /// Function types register a factory lambda under their name(), which is then
        /// invoked to instantiate them. New types are added without touching core code.
        class RegFunctions
        {
        public:
            using func_uptr = std::unique_ptr<Function>;

            /// @brief Check whether a factory is registered under the given name.
            bool contains(std::string_view name) const { return factory_map_.count(name); }

            /// @brief Get all registered names.
            std::vector<std::string> names() const;

            /// @brief Create a Function instance from its registered name.
            /// @param name The name to instantiate.
            /// @param is_prim Whether the object evaluates as its primitive.
            /// @return A unique_ptr to the created object. nullptr if the name is unknown.
            virtual func_uptr create(std::string_view name, bool is_prim = false) const;

            /// @brief Register a function factory under its name.
            /// @tparam T The function class type to register.
            template <typename T>
            void register_factory()
            {
                factory_map_.emplace(std::string(T().name()),
                                     []() -> std::unique_ptr<T>
                                     { return std::make_unique<T>(); });
            }

        protected:
            /// @brief Name to factory map, keyed by the registered names.
            std::map<std::string, std::function<func_uptr()>, std::less<>> factory_map_ = {};
        };

        /// @class RegDistributions
        /// @brief Name to factory map for Distribution objects.
        class RegDistributions : public RegFunctions
        {
        public:
            /// @brief Create a Distribution instance from its registered name.
            /// @param name The name to instantiate.
            /// @param is_cumul Whether the object evaluates as its cumulative.
            /// @return A unique_ptr to the created object. nullptr if the name is unknown.
            func_uptr create(std::string_view name, bool is_cumul = false) const override final;
        };

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

        /// Registry holding the given kind.
        /// NOTE: I define instance in source to avoid duplication across DLLs.
        RegFunctions &instance(Kind kind);

        /// Distribution registries in lookup order. 
        /// NOTE: names shared by both registries resolves to discrete when the distribution is discrete.
        inline std::vector<Kind> distr_kinds(bool discr)
        {
            return discr ? std::vector<Kind>{Kind::Discrete, Kind::Continuous}
                         : std::vector<Kind>{Kind::Continuous, Kind::Discrete};
        }

        /// Create a function from the first registry carrying the name.
        /// An NA name carries no definition and yields a null pointer.
        std::unique_ptr<Function> create_function(
            std::string_view name,
            ErrType err_type,
            const std::vector<Kind> &kinds,
            bool is_prim = false);

        /// Create a distribution from the first kind carrying the name.
        /// An NA name carries no definition and yields a null pointer.
        std::unique_ptr<Distribution> create_distribution(
            std::string_view name,
            ErrType err_type,
            const std::vector<Kind> &kinds,
            bool is_cumul = false);

        /// List the names registered across the given kinds.
        std::vector<std::string> list_names(const std::vector<Kind> &kinds);
    }
}

/// Define the name a function class registers under.
#define HOUSEPI_OVERRIDE_NAME_METHOD(name_) \
    std::string name() const noexcept override final { return #name_; }

/// Define the copy of a function class.
#define HOUSEPI_OVERRIDE_FUNCTION_CLONE_METHOD(Class_) \
    std::unique_ptr<housepi::Function> clone() const override final { return std::make_unique<Class_>(*this); }

/// Define the copy of a distribution class.
#define HOUSEPI_OVERRIDE_DISTRIBUTION_CLONE_METHOD(Class_) \
    std::unique_ptr<housepi::Distribution> clone_distrib() const override final { return std::make_unique<Class_>(*this); }

/// Register a function class into the registry of the given kind, at load time.
#define HOUSEPI_REGISTER_FUNCTION(T, kind)                                 \
    namespace                                                              \
    {                                                                      \
        struct T##kind##_AutoRegister                                      \
        {                                                                  \
            T##kind##_AutoRegister()                                       \
            {                                                              \
                housepi::registry::instance(housepi::registry::Kind::kind) \
                    .register_factory<T>();                                \
            }                                                              \
        };                                                                 \
        static T##kind##_AutoRegister T##kind##_auto_register;             \
    }

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////  *** registry : models ***  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

namespace housepi::registry
{
    /// @class RegModels
    /// @brief Name to factory map for ModelBase objects.
    ///
    /// Model DLLs auto-register their factory into this singleton at load time.
    /// The singleton lives in housepi.so, so all model DLLs share the same instance.
    class RegModels
    {
    public:
        using model_uptr = std::unique_ptr<ModelBase>;

        /// @brief Get the global model registry.
        /// NOTE: I define instance in source to avoid duplication across DLLs.
        static RegModels &instance();

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
        model_uptr create(std::string_view name) const;

    private:
        std::map<std::string, std::function<model_uptr()>, std::less<>> factory_map_;
    };
}

/// Register a module model under the module name, at load time.
#define HOUSEPI_REGISTER_MODEL(module)                                  \
    namespace                                                           \
    {                                                                   \
        struct module##_ModelAutoRegister                               \
        {                                                               \
            module##_ModelAutoRegister()                                \
            {                                                           \
                housepi::registry::RegModels::instance()                \
                    .register_factory<housepi::module::Model>(#module); \
            }                                                           \
        };                                                              \
        static module##_ModelAutoRegister module##_model_auto_register; \
    }
