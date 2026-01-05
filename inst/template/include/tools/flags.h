/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__FLAGS__H
#define HOUSEPI__TOOLS__FLAGS__H

#include <string>

/// Compile-time list of standard housepi flags.
class HousepiRcppCompiledFlags
{
public:
    HousepiRcppCompiledFlags();

    /// Get compilation flags.
    std::string get_flags() { return flags_; }

private:
    std::string flags_;
};

#endif // HOUSEPI__TOOLS__FLAGS__H
