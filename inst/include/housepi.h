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

#include "individual.h"
#include "household.h"

#ifdef _HOUSEPI_RCPP
#include <Rcpp.h>
#else
class MCMC;
int standalone(MCMC &model, int argc, char *argv[]);
#endif

#ifndef _HOUSEPI_MODULE
#define _HOUSEPI_MODULE blank
#endif

#define HOUSEPI_DECLARE_EMPTY_KEYS(T)                        \
    namespace                                                \
    {                                                        \
        static const std::vector<std::string> T##_keys = {}; \
    }

#define HOUSEPI_DECLARE_REGISTER_KEYS(T, ...)                                                        \
    namespace                                                                                        \
    {                                                                                                \
        struct T                                                                                     \
        {                                                                                            \
            enum : size_t                                                                            \
            {                                                                                        \
                __VA_ARGS__                                                                          \
            };                                                                                       \
        };                                                                                           \
                                                                                                     \
        static const std::vector<std::string> T##_keys = hsp::key::helper_parse_names(#__VA_ARGS__); \
    }

#define HOUSEPI_MODEL_CONSTRUCTOR()     \
    Model() : MCMC(DataKey_keys,        \
                   HHsizeKey_keys,      \
                   HHsizevecKey_keys,   \
                   HHdoublevecKey_keys, \
                   HHdoublemtxKey_keys, \
                   LongKey_keys,        \
                   ParamKey_keys,       \
                   GroupRateKey_keys,   \
                   AugmKey_keys)        \
    {                                   \
        init_advanced_setup();          \
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
#define HOUSEPI_RCPP_MODEL_FACTORY()          \
    Rcpp::XPtr<MCMC> housepi_new_model()      \
    {                                         \
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
#define HOUSEPI_REGISTER_RCPP_MODULE_INTERFACE()   \
    int main(int argc, char *argv[])               \
    {                                              \
        using namespace _HOUSEPI_MODULE;           \
        Model user_model;                          \
        return standalone(user_model, argc, argv); \
    }
#endif
