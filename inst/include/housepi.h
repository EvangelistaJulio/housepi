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

#ifdef _STANDALONE
#include "mcmc.h"
#include "standalone.h"
#else
#include <Rcpp.h>
#include "tools/plugin.h"
#endif

#include "model.h"
#include "individual.h"
#include "household.h"

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace hsp::key
{
    /// Parse comma-separated enum names from stringified variadic arguments.
    std::vector<std::string> helper_parse_names(std::string input);
}

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
        return hsp::key::helper_parse_names(#__VA_ARGS__);        \
    }

#define HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE_IMPL_(module_name) \
    namespace hsp::module_name                                 \
    {

#ifdef _HOUSEPI_MODULE
#define HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE(module_name) \
    HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE_IMPL_(_HOUSEPI_MODULE)
#else
#define HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE(module_name) \
    HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE_IMPL_(module_name)
#endif

#define HOUSEPI_CLOSE_NAMED_MODULE_NAMESPACE() }

#ifdef _HOUSEPI_DISCRETE_TIME
#define FLAG_STR_HOUSEPI_DISCRETE_TIME "_HOUSEPI_DISCRETE_TIME,"
#else
#define FLAG_STR_HOUSEPI_DISCRETE_TIME ""
#endif

#ifdef _STANDALONE
#define HOUSEPI_EXTRA_FLAG_SUPPORT()
#define HOUSEPI_RCPP_MODEL_FACTORY()
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(module_name) \
    int main(int argc, char *argv[])                              \
    {                                                             \
        using namespace hsp::module_name;                         \
        Model model = Model();                                    \
        MCMC engine(model);                                       \
        return standalone(engine, argc, argv);                    \
    }
#else
#define HOUSEPI_EXTRA_FLAG_SUPPORT()                                         \
    std::vector<std::string> compile_flags() const override final            \
    {                                                                        \
        return hsp::key::helper_parse_names(FLAG_STR_HOUSEPI_DISCRETE_TIME); \
    }

#define HOUSEPI_RCPP_MODEL_FACTORY()         \
    ModelXPtr new_model()                    \
    {                                        \
        return ModelXPtr(new Model(), true); \
    }

#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(module_name) \
    RCPP_MODULE(module_name)                                      \
    {                                                             \
        using namespace hsp::module_name;                         \
                                                                  \
        Rcpp::function("new_model", &new_model);                  \
    }
#endif

#ifdef _HOUSEPI_MODULE
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE(module_name) \
    HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(_HOUSEPI_MODULE)
#else
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE(module_name) \
    HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(module_name)
#endif
