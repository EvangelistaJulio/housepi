/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <vector>
#include <map>

namespace housepi
{
    template <typename T>
    using dict = std::map<std::string, T, std::less<>>;
    using str_table = dict<std::vector<std::string>>;
    using dbl_table = dict<std::vector<double>>;
}
