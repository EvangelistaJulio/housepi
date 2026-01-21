/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include "tools/keys.h"
#include "tools/keys_helper.h"

#include "individual.h"
#include "household.h"

/// Define enum keys and auto-generate registration.
#define REGISTER_HOUSEPI_KEYS(T, ...)                                    \
    enum class T : size_t                                                \
    {                                                                    \
        __VA_ARGS__                                                      \
    };                                                                   \
                                                                         \
    namespace                                                            \
    {                                                                    \
        struct T##_AutoRegister                                          \
        {                                                                \
            T##_AutoRegister()                                           \
            {                                                            \
                auto keys = hsp::key::helper_parse_names(#__VA_ARGS__);  \
                hsp::key::Registry<T>::instance().init(std::move(keys)); \
            }                                                            \
        };                                                               \
        static T##_AutoRegister T##_housepi_auto_register;               \
    }
