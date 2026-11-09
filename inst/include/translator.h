/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>

#include "tools/dictionary.h"
#include "function.h"

//////////////////////////////  Translator  //////////////////////////////

namespace housepi
{
    namespace parse
    {
        struct Entry;
    }
    namespace rates
    {
        struct Group;
    }
    struct Augmented;
    struct Data;
    struct Parameter;
    class Individual;
    class Household;
    class Variable;

    /// Translate parsed definitions into Attributes, across modules.
    struct Translator
    {
        /// Bind the modules sources are read from.
        /// * `data`, `group_rate`: read by member-specific sources.
        Translator(Parameter &param,
                   const Augmented &augm,
                   const Data &data,
                   const rates::Group &group_rate)
            : param_(param),
              augm_(augm),
              data_(data),
              group_rate_(group_rate) {}

        /// Translate definition node `node` into an attribute.
        /// * within household `house_ptr`, no member bound -> a data cell is read through a fold only.
        const Attribute &to_attribute(const parse::Entry &node,
                                      const Individual *indiv_ptr = nullptr,
                                      const Household *house_ptr = nullptr) const;
        /// Translate expression `expr` into an attribute.
        const Attribute &to_attribute(const std::string &expr,
                                      const Individual *indiv_ptr = nullptr) const;

        /// Translate definition node `node` into a numeric constant attribute.
        const Attribute &to_numeric(const parse::Entry &node,
                                    const Individual *indiv_ptr = nullptr) const;

        /// Translate the arguments of call `node` into numeric constants, keyed by name.
        dict<Attribute> to_numeric_args(const parse::Entry &node,
                                        const Individual *indiv_ptr = nullptr) const;
        /// Bind the arguments of distribution `model` and the bounds `min`, `max` into `var`, then reshape it.
        /// * `numeric` -> every argument a numeric constant;
        /// * otherwise -> arguments live, read on `indiv_ptr`.
        /// @note Call after Variable::init_distr().
        void bind_distr(Variable &var,
                        const parse::Entry &model,
                        const parse::Entry &min,
                        const parse::Entry &max,
                        bool discr,
                        const Individual *indiv_ptr,
                        bool numeric) const;

        /// Mark parameter node `node` as a producer of the transmission cache.
        void mark_param_produces_transm_cache(const parse::Entry &node) const;

        /// Whether reduction `node` folds over the members of household `house_ptr`.
        bool is_house_fold(const parse::Entry &node,
                           const Household *house_ptr) const;
        /// Bind the arguments of function `node` into `parent`.
        /// * household fold -> each operand reading a member cell binds once per member, empty ones skipped.
        /// * defined call: every operand outside a fold carries a value, each fold holds at least one.
        /// @return Whether the call is defined within the household.
        bool bind_arguments(Attribute &parent,
                            const parse::Entry &node,
                            const Individual *indiv_ptr,
                            const Household *house_ptr) const;

    private:
        /// Parameter key module.
        Parameter &param_;
        /// Latent key module.
        const Augmented &augm_;
        /// Data tables and key holder.
        const Data &data_;
        /// Group-rate storage and key holder.
        const rates::Group &group_rate_;
        /// Auxiliary storage attribute.
        mutable Attribute storage_;

        /// Build function `node` into `storage_`, its arguments bound.
        void helper_store_read_function_(const parse::Entry &node,
                                         const Individual *indiv_ptr,
                                         const Household *house_ptr) const;
        /// Whether `node` reads a member value (data cell, latent or group rate) outside an inner reduction.
        bool reads_member_cell_(const parse::Entry &node) const;
        bool reads_member_cell_(const std::vector<parse::Entry> &nodes) const;
    };
}
