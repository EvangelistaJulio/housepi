/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__MESSAGES__H
#define HOUSEPI__TOOLS__MESSAGES__H

#include <string>

namespace hsp
{
    namespace msg
    {
        /// Throw an error message.
        [[noreturn]] void error(const std::string &txt);
        /// Show a message to stdout and flush.
        void info(const std::string &txt);
        /// Display a warning message if not already displayed before.
        void warning(const std::string &txt);
    }
}

#endif // HOUSEPI__TOOLS__MESSAGES__H
