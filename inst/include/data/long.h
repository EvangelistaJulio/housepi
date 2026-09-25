/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

#include "tools/keys.h"
#include "tools/types.h"

#include "data/individual.h"

//////////////////////////////  data::Long  //////////////////////////////

namespace housepi::data
{
    struct Wide;

    /// Long-format table: any number of rows per member.
    struct Long : public KeyModule
    {
        std::vector<std::string> headers = {}; /// User column names, `id` excluded.

        /// Rows of one member.
        struct IndivLong
        {
            std::vector<std::vector<double>> vec = {}; /// Values per column, one per row.
        };
        std::vector<IndivLong> indiv_datalist = {}; /// Member rows, by member index.

        //////  Setter  //////

        /// Build member rows from long-format table `input`, matched on `id` to `data_wide`.
        void build(Wide &data_wide, const str_table &input);

        //////  Getter  //////

        /// Rows of member `indiv`.
        const std::vector<std::vector<double>> *get_data_long_ptr(const Individual &indiv) const noexcept
        {
            return &indiv_datalist[Index::of(indiv)].vec;
        }

    private:
        /// Clear the previous rows.
        void helper_build_clear_data_();
        /// Throw if the headers of `input` are malformed.
        void helper_build_data_headers_(const Wide &data_wide, const str_table &input);
    };
}
