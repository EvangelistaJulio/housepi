/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/dictionary.h"
#include "tools/keys.h"
#include "data/individual.h"

//////////////////////////////  data::Long  //////////////////////////////

namespace housepi::data
{
    struct Wide;

    /// Long-format table: any number of rows per member.
    struct Long : public KeyHolder
    {
        using KeyHolder::KeyHolder;

        /// Rows of one member: values per column, one per row.
        using IndivLong = std::vector<std::vector<double>>;

        std::vector<std::string> headers = {};      /// User column names, `id` excluded.
        std::vector<IndivLong> indiv_datalist = {}; /// Member rows, by member index.

        /// Build member rows from long-format table `input`, matched on `id` to `data_wide`.
        void build(Wide &data_wide, const str_table &input);

        /// Rows of member `indiv`.
        const IndivLong *get_data_long_ptr(const Individual &indiv) const noexcept
        {
            return &indiv_datalist[Index::of(indiv)];
        }

    private:
        /// Clear the previous rows.
        void helper_build_clear_data_();
        /// Throw if the headers of `input` are malformed.
        void helper_build_data_headers_(const Wide &data_wide,
                                        const str_table &input);
    };
}
