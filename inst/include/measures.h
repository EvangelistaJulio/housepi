/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>

namespace housepi::measure
{
    namespace cases
    {
        struct Report
        {
            struct Summary
            {
                unsigned house_size = 0;
                unsigned cases = 0;
            };

            std::vector<Summary> houses;
            unsigned min_size = 0;
            unsigned max_size = 0;
        };

        struct Table
        {
            unsigned min_size = 0;
            unsigned max_size = 0;
            /// counts[row][c] with row = house_size - min_size, columns 0..max_size.
            std::vector<std::vector<unsigned>> counts;
        };

        struct NonZero
        {
            unsigned max_size = 0;       /// 0 if no qualifying households.
            std::vector<unsigned> count; /// Sizes 2..max_size; empty if max_size < 2.
        };
    }

    struct SecAttackRates
    {
        unsigned grid_min = 0;        /// First house size on the grid (>= 2).
        std::vector<double> sar_mean; /// Mean SAR per size; 0.0 where count == 0.
        std::vector<unsigned> count;  /// Number of contributing households per size.
    };

    /// Aggregate individual statuses by house_id, tracking the size range in the same pass.
    /// `status` values: negative = missing/NA, 0 = susceptible, >=1 = infected.
    cases::Report summarise_houses(
        const std::vector<std::string> &house_id,
        const std::vector<int> &status);

    SecAttackRates sar_by_house_size(const cases::Report &report);
    cases::Table cases_by_house_size(const cases::Report &report);
    cases::NonZero non_zero_cases(const cases::Report &report);
}
