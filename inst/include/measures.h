/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <algorithm>
#include <vector>

#include "data/status.h"

namespace housepi::measure
{
    namespace cases
    {
        /// Count households per size and number of cases.
        struct Table
        {
            /// Smallest household size observed.
            unsigned min_size = 0;
            /// Largest household size observed.
            unsigned max_size = 0;
            /// Rows: hh_size - min_size; columns: cases 0..max_size.
            std::vector<std::vector<unsigned>> counts;
        };

        /// Count households per size, with at least one case or in total.
        struct HouseCounts
        {
            /// First household size on the grid, 1 or more.
            unsigned grid_min = 0;
            /// Households per size.
            std::vector<unsigned> total;
            /// Households with at least one case, per size.
            std::vector<unsigned> at_least_1_case;
        };

        /// Hold mean attack rate per household size.
        struct AttackRates
        {
            /// First household size on the grid, 1 or more.
            unsigned grid_min = 0;
            /// Mean attack rate per size.
            std::vector<double> ar_mean;
            /// Contributing households per size.
            std::vector<unsigned> count;
        };

        /// Hold mean secondary attack rate per household size.
        struct SecAttackRates
        {
            /// First household size on the grid, 2 or more.
            unsigned grid_min = 0;
            /// Mean secondary attack rate per size.
            std::vector<double> sar_mean;
            /// Contributing households per size.
            std::vector<unsigned> count;
        };
    }

    /// Tabulate households by size and number of cases, aggregating individuals by `hh_id`.
    ///
    /// * `hh_id`: 1-based codes 1..n_houses, each code used at least once;
    /// * `status`: values as read from an R column, NA missing, 0 escaped, 1 or more infected.
    inline cases::Table summarise_houses(const std::vector<int> &hh_id,
                                         const std::vector<int> &status)
    {
        const size_t n_indiv = status.size();
        if (n_indiv == 0)
            return cases::Table();

        struct House
        {
            unsigned hh_size = 0;
            unsigned cases = 0;
        };
        std::vector<House> houses(*std::max_element(hh_id.begin(), hh_id.end()));
        for (size_t i = 0; i < n_indiv; ++i)
        {
            House &h = houses[hh_id[i] - 1];
            ++h.hh_size;
            h.cases += status::is_infected(status[i]);
        }

        cases::Table output;
        output.min_size = -1;
        for (const auto &h : houses)
        {
            output.min_size = std::min(output.min_size, h.hh_size);
            output.max_size = std::max(output.max_size, h.hh_size);
        }

        output.counts.assign(output.max_size - output.min_size + 1,
                             std::vector<unsigned>(output.max_size + 1, 0));
        for (const auto &h : houses)
            ++output.counts[h.hh_size - output.min_size][h.cases];

        return output;
    }

    /// Count households per size, all and with a case, on the `attack_rate_by_hh_size` grid.
    inline cases::HouseCounts house_counts(const cases::Table &table)
    {
        cases::HouseCounts output;
        output.grid_min = std::max(table.min_size, 1u);
        for (const auto &row : table.counts)
        {
            unsigned total = 0;
            for (unsigned n_houses : row)
                total += n_houses;
            output.total.push_back(total);
            output.at_least_1_case.push_back(total - row[0]);
        }
        return output;
    }

    /// Compute mean attack rate per household size.
    inline cases::AttackRates attack_rate_by_hh_size(const cases::Table &table)
    {
        cases::AttackRates output;
        output.grid_min = std::max(table.min_size, 1u);
        for (unsigned r = 0; r < table.counts.size(); ++r)
        {
            unsigned total = 0, n_cases = 0;
            for (unsigned c = 0; c < table.counts[r].size(); ++c)
            {
                total += table.counts[r][c];
                n_cases += c * table.counts[r][c];
            }
            output.count.push_back(total);
            output.ar_mean.push_back(
                total > 0 ? static_cast<double>(n_cases) / ((table.min_size + r) * total) : 0.0);
        }
        return output;
    }

    /// Compute mean secondary attack rate per household size, from size 2.
    inline cases::SecAttackRates sar_by_hh_size(const cases::Table &table)
    {
        cases::SecAttackRates output;
        output.grid_min = std::max(table.min_size, 2u);
        for (unsigned r = output.grid_min - table.min_size; r < table.counts.size(); ++r)
        {
            unsigned affected = 0, secondary = 0;
            for (unsigned c = 1; c < table.counts[r].size(); ++c)
            {
                affected += table.counts[r][c];
                secondary += (c - 1) * table.counts[r][c];
            }
            output.count.push_back(affected);
            output.sar_mean.push_back(
                affected > 0 ? static_cast<double>(secondary) / ((table.min_size + r - 1.0) * affected) : 0.0);
        }
        return output;
    }
}
