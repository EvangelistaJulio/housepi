/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__SOURCES__H
#define HOUSEPI__TOOLS__SOURCES__H

#include <string>

namespace hsp
{
    // Source checking functions.
    namespace sources
    {
        inline bool is_augmented(const std::string &s) { return s == "A"; }
        inline bool is_data(const std::string &s) { return s == "D"; }
        inline bool is_group_rate(const std::string &s) { return s == "G"; }
        inline bool is_numeric(const std::string &s) { return s == "N"; }
        inline bool is_parameter(const std::string &s) { return s == "P"; }
    }
}

#endif // HOUSEPI__TOOLS__SOURCES__H
