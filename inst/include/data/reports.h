/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstddef>
#include <utility>

//////////////////////////////  data::Report*  //////////////////////////////

namespace housepi::data
{
    /// Counts over every member.
    struct ReportIndividual
    {
        size_t total = 0;      /// Members.
        size_t missing = 0;    /// Members with missing status.
        size_t escaped = 0;    /// Members with escaped status.
        size_t infectious = 0; /// Members with infected status.
    };
    /// Counts over every household.
    struct ReportHousehold
    {
        size_t total = 0;            /// Households.
        size_t distinct = 0;         /// Equivalence classes.
        size_t no_cases = 0;         /// Households with no infected member.
        size_t at_least_1_case = 0;  /// Households with at least one infected member.
        size_t at_least_2_cases = 0; /// Households with at least two infected members.
    };
    /// Extremes over household sizes, each paired with the number of households reaching it.
    struct ReportHouseSizes
    {
        std::pair<unsigned short, size_t> max_cases = {0, 0};   /// Most infected members.
        std::pair<unsigned short, size_t> max_miss = {0, 0};    /// Most members with missing status.
        std::pair<unsigned short, size_t> min_size = {1000, 0}; /// Fewest members.
        std::pair<unsigned short, size_t> max_size = {0, 0};    /// Most members.
    };
}
