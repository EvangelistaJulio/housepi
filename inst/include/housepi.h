/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <vector>

#include "tools/math.h"
#include "tools/messages.h"
#include "tools/string.h"
#include "data/household.h"
#include "data/individual.h"
#include "model.h"
#include "registry.h"

#ifdef _STANDALONE
#include "mcmc.h"
#include "standalone.h"
#endif

/// Declare keys `...` as `T::key`, and override `method` to list their names.
#define HOUSEPI_DECLARE_KEYS_IMPL(T, method, ...)          \
    struct T                                               \
    {                                                      \
        enum : size_t                                      \
        {                                                  \
            __VA_ARGS__                                    \
        };                                                 \
    };                                                     \
                                                           \
    std::vector<std::string> method() const override final \
    {                                                      \
        return housepi::str::parse_names(#__VA_ARGS__);    \
    }

/// Declare the keys a model reads, one macro per slot, none for a slot it leaves empty.
#define HOUSEPI_DECLARE_DATA_KEYS(...) HOUSEPI_DECLARE_KEYS_IMPL(Data, data_keys, __VA_ARGS__)
#define HOUSEPI_DECLARE_LONG_KEYS(...) HOUSEPI_DECLARE_KEYS_IMPL(LongDataset, long_keys, __VA_ARGS__)
#define HOUSEPI_DECLARE_PARAM_KEYS(...) HOUSEPI_DECLARE_KEYS_IMPL(Param, param_keys, __VA_ARGS__)
#define HOUSEPI_DECLARE_GROUP_RATE_KEYS(...) HOUSEPI_DECLARE_KEYS_IMPL(GroupRate, group_rate_keys, __VA_ARGS__)
#define HOUSEPI_DECLARE_LATENT_KEYS(...) HOUSEPI_DECLARE_KEYS_IMPL(Latent, latent_keys, __VA_ARGS__)
#define HOUSEPI_DECLARE_HOUSE_MEMBER_KEYS(...) \
    HOUSEPI_DECLARE_KEYS_IMPL(HouseMember, house_member_keys, __VA_ARGS__)
#define HOUSEPI_DECLARE_HOUSE_MEMBER_VEC_KEYS(...) \
    HOUSEPI_DECLARE_KEYS_IMPL(HouseMemberVec, house_member_vec_keys, __VA_ARGS__)
#define HOUSEPI_DECLARE_HOUSE_DOUBLE_VEC_KEYS(...) \
    HOUSEPI_DECLARE_KEYS_IMPL(HouseDoubleVec, house_double_vec_keys, __VA_ARGS__)
#define HOUSEPI_DECLARE_HOUSE_DOUBLE_MTX_KEYS(...) \
    HOUSEPI_DECLARE_KEYS_IMPL(HouseDoubleMtx, house_double_mtx_keys, __VA_ARGS__)

#ifdef _STANDALONE
#define HOUSEPI_STANDALONE_MAIN(module)                  \
    int main(int argc, char *argv[])                     \
    {                                                    \
        return housepi::standalone(#module, argc, argv); \
    }
#define HOUSEPI_REGISTER_MODULE_INTERFACE_IMPL_(module) \
    HOUSEPI_REGISTER_MODEL(module)                      \
    HOUSEPI_STANDALONE_MAIN(module)

#else
#define HOUSEPI_REGISTER_MODULE_INTERFACE_IMPL_(module) \
    HOUSEPI_REGISTER_MODEL(module)
#endif

#define HOUSEPI_REGISTER_MODULE_INTERFACE() \
    HOUSEPI_REGISTER_MODULE_INTERFACE_IMPL_(_HOUSEPI_MODULE)

#define HOUSEPI_OPEN_MODULE_NAMESPACE_IMPL_(module)    \
    namespace housepi::HOUSEPI_MODEL_NAMESPACE(module) \
    {
#define HOUSEPI_OPEN_MODULE_NAMESPACE() \
    HOUSEPI_OPEN_MODULE_NAMESPACE_IMPL_(_HOUSEPI_MODULE)
#define HOUSEPI_CLOSE_MODULE_NAMESPACE() \
    }                                    \
    HOUSEPI_REGISTER_MODULE_INTERFACE()
