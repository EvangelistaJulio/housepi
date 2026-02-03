/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include "tools/keys.h"

#include "individual.h"
#include "household.h"

#define EMPTY_HOUSEPI_KEYS(T)                                \
    namespace                                                \
    {                                                        \
        static const std::vector<std::string> T##_keys = {}; \
    }

#define REGISTER_HOUSEPI_KEYS(T, ...)                                                                \
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

#define HOUSEPI_MODEL_CONSTRUCTOR()            \
    Model::Model() : MCMC(DataKey_keys,        \
                          HHsizeKey_keys,      \
                          HHsizevecKey_keys,   \
                          HHdoublevecKey_keys, \
                          HHdoublemtxKey_keys, \
                          LongKey_keys,        \
                          ParamKey_keys,       \
                          GroupRateKey_keys,   \
                          AugmKey_keys)        \
    {                                          \
        init_advanced_setup();                 \
    }
