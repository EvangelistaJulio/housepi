/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "tools/status.h"

namespace housepi::measure
{
    namespace cases
    {
        struct Table
        {
            /// Smallest household size observed.
            unsigned min_size = 0;
            /// Largest household size observed.
            unsigned max_size = 0;
            /// Most cases observed in a single household.
            unsigned max_cases = 0;
            /// counts[r][c] with rows = hhsize - min_size, columns 0..max_size.
            std::vector<std::vector<unsigned>> counts;
        };

        struct HouseCounts
        {
            /// First house size on the grid (>= 1).
            unsigned grid_min = 0;
            /// All houses per size from grid_min.
            std::vector<unsigned> total;
            /// Houses with >=1 case per size.
            std::vector<unsigned> at_least_1_case;
        };

        struct AttackRates
        {
            /// First house size on the grid (>= 1).
            unsigned grid_min = 0;
            /// Mean attack rate per size.
            std::vector<double> ar_mean;
            /// Contributing households per size.
            std::vector<unsigned> count;
        };

        struct SecAttackRates
        {
            /// First house size on the grid (>= 2).
            unsigned grid_min = 0;
            /// Mean SAR per size.
            std::vector<double> sar_mean;
            /// Contributing households per size.
            std::vector<unsigned> count;
        };
    }

    /// Tabulate households by size and number of cases, aggregating individuals by house_id.
    /// `house_id` holds 1-based codes 1..n_houses, each code used at least once.
    /// `status` values as read from an R column: NA = missing, 0 = escaped, >=1 = infected.
    inline cases::Table summarise_houses(const std::vector<int> &house_id,
                                         const std::vector<int> &status)
    {
        const size_t n_indiv = status.size();
        if (n_indiv == 0)
            return cases::Table();

        struct House
        {
            unsigned hhsize = 0;
            unsigned cases = 0;
        };
        std::vector<House> houses(*std::max_element(house_id.begin(), house_id.end()));
        for (size_t i = 0; i < n_indiv; ++i)
        {
            House &h = houses[house_id[i] - 1];
            ++h.hhsize;
            h.cases += status::is_infected(status[i]);
        }

        cases::Table output;
        output.min_size = -1;
        for (const auto &h : houses)
        {
            output.min_size = std::min(output.min_size, h.hhsize);
            output.max_size = std::max(output.max_size, h.hhsize);
            output.max_cases = std::max(output.max_cases, h.cases);
        }

        output.counts.assign(output.max_size - output.min_size + 1,
                             std::vector<unsigned>(output.max_size + 1, 0));
        for (const auto &h : houses)
            ++output.counts[h.hhsize - output.min_size][h.cases];

        return output;
    }

    /// Total houses and houses with >=1 case per size, on the same grid as attack_rate_by_hhsize.
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

    inline cases::AttackRates attack_rate_by_hhsize(const cases::Table &table)
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

    inline cases::SecAttackRates sar_by_hhsize(const cases::Table &table)
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

    /// Geometric mean of household size, weighted by secondary cases max(cases - 1, 0).
    /// NOTE: requires cases::Table::max_cases > 1.
    inline double geom_mean_hhsize(const cases::Table &table)
    {
        double sum_w = 0.0, sum_wlog = 0.0;
        for (unsigned r = 0; r < table.counts.size(); ++r)
        {
            unsigned w = 0;
            for (unsigned c = 2; c < table.counts[r].size(); ++c)
                w += (c - 1) * table.counts[r][c];
            sum_w += w;
            sum_wlog += w * std::log(static_cast<double>(table.min_size + r));
        }
        return std::exp(sum_wlog / sum_w);
    }
}
