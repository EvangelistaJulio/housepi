/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>

#include "tools/messages.h"
#include "model.h"

//////////////////////////////  registry::RegModels  //////////////////////////////

namespace housepi::registry
{
    /// Map names to factories of ModelBase objects.
    ///
    /// * model DLLs register their factory into this singleton at load time;
    /// * a name registered again -> the last loaded factory;
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
            factory_map_.insert_or_assign(
                name,
                Factory{[]() -> model_uptr
                        { return std::make_unique<T>(); },
                        [](const ModelBase &src) -> model_uptr
                        { return std::make_unique<T>(static_cast<const T &>(src)); }});
        }

        /// Create the model registered under `name`.
        /// * unregistered name -> error.
        model_uptr create(std::string_view name) const { return find_(name).create(); }
        /// Copy `src`, a model created under `name`.
        /// * unregistered name -> error.
        model_uptr copy(std::string_view name, const ModelBase &src) const { return find_(name).copy(src); }

    private:
        /// Factories of one model class.
        struct Factory
        {
            /// New model.
            std::function<model_uptr()> create;
            /// Copy of a model of the class.
            std::function<model_uptr(const ModelBase &)> copy;
        };

        /// Factories, keyed by model name.
        std::map<std::string, Factory, std::less<>> factory_map_;

        /// Factories registered under `name`.
        /// * unregistered name -> error.
        const Factory &find_(std::string_view name) const
        {
            auto it = factory_map_.find(name);
            // # nocov start
            if (it == factory_map_.end())
                msg::stop("Model registration failed: unknown model " + msg::val(std::string(name)) + ".\n");
            // # nocov end
            return it->second;
        }
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
