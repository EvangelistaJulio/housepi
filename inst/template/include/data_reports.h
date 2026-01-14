/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__DATA_REPORTS__H
#define HOUSEPI__DATA_REPORTS__H

#include <utility>
#include <vector>

/// Report for all individuals in the data.
struct ReportIndividual
{
	size_t total = 0;
	size_t missing = 0;
	size_t susceptible = 0;
	size_t infectious = 0;
};
/// Report for all households in the data.
struct ReportHousehold
{
	size_t total = 0;
	size_t no_cases = 0;
	size_t at_least_1_case = 0;
	size_t at_least_2_cases = 0;
};
/// Report for all household sizes.
struct ReportHouseSizes
{
	std::pair<unsigned short, size_t> max_cases = {0, 0};
	std::pair<unsigned short, size_t> max_miss = {0, 0};
	std::pair<unsigned short, size_t> min_size = {1000, 0};
	std::pair<unsigned short, size_t> max_size = {0, 0};
};

#endif // HOUSEPI__DATA_REPORTS__H
