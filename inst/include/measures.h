/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "tools/status.h"

namespace housepi::measure
{
    namespace cases
    {
        struct Report
        {
            struct Summary
            {
                unsigned hhsize = 0;
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
            /// counts[row][c] with row = hhsize - min_size, columns 0..max_size.
            std::vector<std::vector<unsigned>> counts;
        };

        struct HouseCounts
        {
            unsigned grid_min = 0;                 /// First house size on the grid (>= 1).
            std::vector<unsigned> total;           /// All houses per size from grid_min.
            std::vector<unsigned> at_least_1_case; /// Houses with >=1 case per size.
        };
    }

    struct AttackRates
    {
        unsigned grid_min = 0;       /// First house size on the grid (>= 1).
        std::vector<double> ar_mean; /// Mean attack rate per size.
        std::vector<unsigned> count; /// N. of contributing households per size.
    };

    struct SecAttackRates
    {
        unsigned grid_min = 0;        /// First house size on the grid (>= 2).
        std::vector<double> sar_mean; /// Mean SAR per size.
        std::vector<unsigned> count;  /// N. of contributing households per size.
    };

    inline cases::Report helper_summarise_single_individuals_(
        const std::vector<int> &status)
    {
        cases::Report output;
        output.houses.reserve(status.size());

        for (int s : status)
            output.houses.push_back({1u, status::is_infected(s)});

        output.min_size = 1;
        output.max_size = 1;

        return output;
    }

    /// Aggregate individual statuses by house_id, tracking the size range in the same pass.
    /// `status` values: negative = missing/NA, 0 = escaped, >=1 = infected.
    inline cases::Report summarise_houses(
        const std::vector<std::string> &house_id,
        const std::vector<int> &status)
    {
        const size_t n_indiv = status.size();

        if (n_indiv == 0)
            return cases::Report();

        if (house_id.size() == 0)
            return helper_summarise_single_individuals_(status);

        cases::Report output;
        output.houses.reserve(n_indiv);

        std::unordered_map<std::string, size_t> id_map;
        for (size_t i = 0; i < n_indiv; ++i)
        {
            auto it = id_map.find(house_id[i]);
            if (it == id_map.end())
            {
                id_map.emplace(house_id[i], output.houses.size());
                output.houses.push_back({1u, status::is_infected(status[i])});
            }
            else
            {
                ++output.houses[it->second].hhsize;
                output.houses[it->second].cases += status::is_infected(status[i]);
            }
        }

        output.min_size = -1;
        for (const auto &h : output.houses)
            if (h.hhsize < output.min_size)
                output.min_size = h.hhsize;
            else if (h.hhsize > output.max_size)
                output.max_size = h.hhsize;

        return output;
    }

    /// Total houses and houses with >=1 case per size, on the same grid as
    /// attack_rate_by_hhsize (grid_min = max(min_size, 1)).
    inline cases::HouseCounts house_counts(const cases::Report &report)
    {
        unsigned grid_min = std::max(report.min_size, 1u);
        unsigned grid_len = report.max_size < grid_min ? 0
                                                       : report.max_size - grid_min + 1;

        cases::HouseCounts output;
        output.grid_min = grid_min;
        if (grid_len == 0)
            return output;

        output.total.assign(grid_len, 0);
        output.at_least_1_case.assign(grid_len, 0);

        for (const auto &h : report.houses)
            if (h.hhsize >= grid_min)
            {
                unsigned idx = h.hhsize - grid_min;
                output.total[idx]++;
                if (h.cases > 0)
                    output.at_least_1_case[idx]++;
            }

        return output;
    }

    inline AttackRates attack_rate_by_hhsize(const cases::Report &report)
    {
        unsigned grid_min = std::max(report.min_size, 1u);
        unsigned grid_len = report.max_size < grid_min ? 0
                                                       : report.max_size - grid_min + 1;

        AttackRates output;
        output.grid_min = grid_min;
        if (grid_len == 0)
            return output;

        output.ar_mean.assign(grid_len, 0.0);
        output.count.assign(grid_len, 0);

        for (const auto &h : report.houses)
            if (h.hhsize >= 1)
            {
                unsigned idx = h.hhsize - grid_min;
                output.ar_mean[idx] += static_cast<double>(h.cases) / h.hhsize;
                output.count[idx]++;
            }
        for (unsigned i = 0; i < grid_len; ++i)
            if (output.count[i] > 0)
                output.ar_mean[i] /= output.count[i];

        return output;
    }

    inline SecAttackRates sar_by_hhsize(const cases::Report &report)
    {
        unsigned grid_min = std::max(report.min_size, 2u);
        unsigned grid_len = report.max_size < grid_min ? 0
                                                       : report.max_size - grid_min + 1;

        SecAttackRates output;
        output.grid_min = grid_min;
        if (grid_len == 0)
            return output;

        output.sar_mean.assign(grid_len, 0.0);
        output.count.assign(grid_len, 0);

        for (const auto &h : report.houses)
            if (h.cases > 0 && h.hhsize > 1)
            {
                unsigned idx = h.hhsize - grid_min;
                output.sar_mean[idx] += static_cast<double>(h.cases - 1.0) / (h.hhsize - 1.0);
                output.count[idx]++;
            }
        for (unsigned i = 0; i < grid_len; ++i)
            if (output.count[i] > 0)
                output.sar_mean[i] /= output.count[i];

        return output;
    }

    inline cases::Table cases_by_hhsize(const cases::Report &report)
    {
        cases::Table output;
        output.min_size = report.min_size;
        output.max_size = report.max_size;

        if (report.houses.empty())
            return output;

        unsigned n_rows = output.max_size - output.min_size + 1;
        output.counts.assign(n_rows, std::vector<unsigned>(output.max_size + 1, 0));
        for (const auto &h : report.houses)
            output.counts[h.hhsize - output.min_size][h.cases]++;

        return output;
    }
}
