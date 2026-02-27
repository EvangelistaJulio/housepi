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

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#include "tools/plugin.h"
#else
#include "mcmc.h"
#include "standalone.h"
#endif

#include "model.h"
#include "individual.h"
#include "household.h"

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE blank
#endif

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

#define HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE() \
    namespace _HOUSEPI_MODULE                 \
    {
#define HOUSEPI_CLOSE_NAMED_MODULE_NAMESPACE() }

/// Expose compile flags.
std::string compile_flags()
{
    return
#ifdef _HOUSEPI_DISCRETE_TIME
        "_HOUSEPI_DISCRETE_TIME:";
#else
        "";
#endif
}

#ifdef _HOUSEPI_RCPP
#define HOUSEPI_RCPP_MODEL_FACTORY()         \
    ModelXPtr new_model()                    \
    {                                        \
        return ModelXPtr(new Model(), true); \
    }

#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(module_name) \
    RCPP_MODULE(module_name)                                      \
    {                                                             \
        using namespace module_name;                              \
                                                                  \
        Rcpp::function("compile_flags", &compile_flags);          \
        Rcpp::function("new_model", &new_model);                  \
    }

#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE() \
    HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(_HOUSEPI_MODULE)
#else
#define HOUSEPI_RCPP_MODEL_FACTORY()
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE() \
    int main(int argc, char *argv[])             \
    {                                            \
        using namespace _HOUSEPI_MODULE;         \
        Model model = Model();                   \
        MCMC engine(model);                      \
        return standalone(engine, argc, argv);   \
    }
#endif
