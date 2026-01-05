/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include "tools/flags.h"

/// Compile-time list of standard housepi flags.
HousepiRcppCompiledFlags::HousepiRcppCompiledFlags()
{
    flags_ = 
#ifdef _HOUSEPI_DISCRETE_TIME
        "_HOUSEPI_DISCRETE_TIME:"
#endif
    "";
}
