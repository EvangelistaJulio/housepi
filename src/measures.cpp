/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include "measures.h"

#include <algorithm>
#include <unordered_map>

#include "tools/status.h"

namespace hsp::measure
{
    cases::Report helper_summarise_single_individuals_(
        const std::vector<int> &status)
    {
        cases::Report output;
        output.houses.reserve(status.size());

        for (int s : status)
            output.houses.push_back({1u, status::is_infected(s)});

        output.min_size = 1;
        output.max_size = 1;

        return output;
    };

    cases::Report summarise_houses(
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
                ++output.houses[it->second].house_size;
                output.houses[it->second].cases += status::is_infected(status[i]);
            }
        }

        output.min_size = -1;
        for (const auto &h : output.houses)
            if (h.house_size < output.min_size)
                output.min_size = h.house_size;
            else if (h.house_size > output.max_size)
                output.max_size = h.house_size;

        return output;
    }

    SecAttackRates sar_by_house_size(const cases::Report &report)
    {
        unsigned grid_min = std::max(report.min_size, 2u);
        unsigned grid_len = report.max_size >= grid_min ? report.max_size - grid_min + 1 : 0;

        SecAttackRates output;
        output.grid_min = grid_min;
        if (grid_len == 0)
            return output;

        output.sar_mean.assign(grid_len, 0.0);
        output.count.assign(grid_len, 0);

        for (const auto &h : report.houses)
            if (h.cases > 0 && h.house_size > 1)
            {
                unsigned idx = h.house_size - grid_min;
                output.sar_mean[idx] += static_cast<double>(h.cases - 1.0) / (h.house_size - 1.0);
                output.count[idx]++;
            }
        for (unsigned i = 0; i < grid_len; ++i)
            if (output.count[i] > 0)
                output.sar_mean[i] /= output.count[i];

        return output;
    }

    cases::Table cases_by_house_size(const cases::Report &report)
    {
        cases::Table output;
        output.min_size = report.min_size;
        output.max_size = report.max_size;

        if (report.houses.empty())
            return output;

        unsigned n_rows = output.max_size - output.min_size + 1;
        output.counts.assign(n_rows, std::vector<unsigned>(output.max_size + 1, 0));
        for (const auto &h : report.houses)
            output.counts[h.house_size - output.min_size][h.cases]++;

        return output;
    }

    cases::NonZero non_zero_cases(const cases::Report &report)
    {
        cases::NonZero output;
        output.max_size = report.max_size;

        if (output.max_size < 2)
            return output;

        output.count.assign(output.max_size - 1, 0); // sizes 2 to max_size
        for (const auto &h : report.houses)
            if (h.cases > 0 && h.house_size > 1)
                output.count[h.house_size - 2]++;

        return output;
    }
}
