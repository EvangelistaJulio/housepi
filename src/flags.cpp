/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#include <string>

std::string housepi_compile_flags()
{
    return
#ifdef _HOUSEPI_DISCRETE_TIME
        "_HOUSEPI_DISCRETE_TIME:";
#else
        "";
#endif
}
