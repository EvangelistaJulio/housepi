/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>

namespace hsp::msg
{
    /// Throw an error message.
    [[noreturn]] void error(const std::string &txt);
    /// Show a message to stdout and flush.
    void info(const std::string &txt) noexcept;
    /// Display a warning message if not already displayed before.
    void warning(const std::string &txt);
    /// Wheter text represent a previously unseen warning.
    bool is_new_housepi_warning(const std::string &txt) noexcept;
    /// Mark given text as an already seen warning.
    void add_new_housepi_warning(const std::string &txt) noexcept;
}
