/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "data/individual.h"

//////////////////////////////  data::LongDataset  //////////////////////////////

namespace housepi::data
{
    struct Dataset;

    /// Long dataset: any number of rows per member.
    struct LongDataset : public KeyHolder
    {
        using KeyHolder::KeyHolder;

        /// Rows of one member: values per column, one per row.
        using IndivRows = std::vector<std::vector<double>>;
        /// Member rows, by member index.
        std::vector<IndivRows> indiv_datalist = {};

        /// Build member rows from data columns `cells`, matched on `ids` column `id` to `dataset`.
        void build(Dataset &dataset, const str_table &ids, const dbl_table &cells);

        /// Rows of member `indiv`.
        const IndivRows *rows_of(const Individual &indiv) const noexcept
        {
            return &indiv_datalist[Index::of(indiv)];
        }

    private:
        /// Clear the previous rows.
        void helper_build_clear_data_();
        /// Throw if `ids` misses column `id`, register the headers of `cells`.
        void helper_build_data_headers_(const Dataset &dataset,
                                        const str_table &ids,
                                        const dbl_table &cells);
    };
}
