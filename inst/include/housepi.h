/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "tools/messages.h"
#include "tools/read.h"
#ifdef _STANDALONE
#include "mcmc.h"
#include "standalone.h"
#endif

#include "model.h"
#include "individual.h"
#include "household.h"
#include "registry.h"

#ifdef _HOUSEPI_DISCRETE_TIME
#define FLAG_STR_HOUSEPI_DISCRETE_TIME "_HOUSEPI_DISCRETE_TIME,"
#else
#define FLAG_STR_HOUSEPI_DISCRETE_TIME ""
#endif

#define HOUSEPI_DECLARE_EMPTY_KEYS(T) \
    struct T                          \
    {                                 \
    };                                \
                                      \
    const std::vector<std::string> get_##T() const override final { return {}; }

#define HOUSEPI_DECLARE_REGISTER_KEYS(T, ...)                     \
    struct T                                                      \
    {                                                             \
        enum : size_t                                             \
        {                                                         \
            __VA_ARGS__                                           \
        };                                                        \
    };                                                            \
                                                                  \
    const std::vector<std::string> get_##T() const override final \
    {                                                             \
        return housepi::read::helper_parse_names(#__VA_ARGS__);       \
    }

#ifdef _STANDALONE
#define HOUSEPI_STANDALONE_MAIN(module)              \
    int main(int argc, char *argv[])                 \
    {                                                \
        return hsp::standalone(#module, argc, argv); \
    }
#define HOUSEPI_EXTRA_FLAG_SUPPORT()
#define HOUSEPI_REGISTER_MODULE_INTERFACE_IMPL_(module) \
    HOUSEPI_REGISTER_MODEL(module)                      \
    HOUSEPI_STANDALONE_MAIN(module)

#else
#define HOUSEPI_EXTRA_FLAG_SUPPORT()                                          \
    std::vector<std::string> compile_flags() const override final             \
    {                                                                         \
        return housepi::read::helper_parse_names(FLAG_STR_HOUSEPI_DISCRETE_TIME); \
    }
#define HOUSEPI_REGISTER_MODULE_INTERFACE_IMPL_(module) \
    HOUSEPI_REGISTER_MODEL(module)
#endif

#define HOUSEPI_REGISTER_MODULE_INTERFACE() \
    HOUSEPI_REGISTER_MODULE_INTERFACE_IMPL_(_HOUSEPI_MODULE)

#define HOUSEPI_OPEN_MODULE_NAMESPACE_IMPL_(module) \
    namespace housepi::module                           \
    {
#define HOUSEPI_OPEN_MODULE_NAMESPACE() \
    HOUSEPI_OPEN_MODULE_NAMESPACE_IMPL_(_HOUSEPI_MODULE)
#define HOUSEPI_CLOSE_MODULE_NAMESPACE() \
    }                                    \
    HOUSEPI_REGISTER_MODULE_INTERFACE()
