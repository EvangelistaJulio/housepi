/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstddef>

//////////////////////////////  data::Report*  //////////////////////////////

namespace housepi::data
{
    /// Counts over every member.
    struct ReportIndividual
    {
        /// Members.
        size_t total = 0;
        /// Members with missing status.
        size_t missing = 0;
        /// Members with escaped status.
        size_t escaped = 0;
        /// Members with infected status.
        size_t infectious = 0;
    };
    /// Counts over every household.
    struct ReportHousehold
    {
        /// Households.
        size_t total = 0;
        /// Households with no infected member.
        size_t no_cases = 0;
        /// Households with at least one infected member.
        size_t at_least_1_case = 0;
        /// Households with at least two infected members.
        size_t at_least_2_cases = 0;
    };
}
