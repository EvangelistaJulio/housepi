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

#include "tools/keys.h"

#include "model.h"
#include "individual.h"
#include "household.h"

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#else
#include "standalone.h"
#endif

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE blank
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
        return hsp::key::helper_parse_names(#__VA_ARGS__);        \
    }

#define HOUSEPI_OPEN_NAMED_MODULE_NAMESPACE() \
    namespace _HOUSEPI_MODULE                 \
    {
#define HOUSEPI_CLOSE_NAMED_MODULE_NAMESPACE() }

/// Expose compile flags.
std::string housepi_compile_flags()
{
    return
#ifdef _HOUSEPI_DISCRETE_TIME
        "_HOUSEPI_DISCRETE_TIME:";
#else
        "";
#endif
}

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#ifdef _HOUSEPI_RCPP
#define HOUSEPI_RCPP_MODEL_FACTORY()                \
    Rcpp::XPtr<MCMC> housepi_new_model()            \
    {                                               \
        return Rcpp::XPtr<MCMC>(new Model(), true); \
    }

#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(module_name)        \
    RCPP_MODULE(module_name)                                             \
    {                                                                    \
        using namespace module_name;                                     \
                                                                         \
        Rcpp::function("housepi_compile_flags", &housepi_compile_flags); \
        Rcpp::function("housepi_new_model", &housepi_new_model);         \
    }

#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE() \
    HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE_IMPL_(_HOUSEPI_MODULE)
#else
#define HOUSEPI_RCPP_MODEL_FACTORY()
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE() \
    int main(int argc, char *argv[])             \
    {                                            \
        using namespace _HOUSEPI_MODULE;         \
        return standalone<Model>(argc, argv);    \
    }
#endif
