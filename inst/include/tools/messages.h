/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <exception>
#include <string>
#include <utility>

namespace housepi
{
    /// Package-authored exception.
    /// NOTE: must keep the name in sync with the R side.
    struct CppError : std::exception
    {
        std::string m;
        explicit CppError(std::string s) : m(std::move(s)) {}
        const char *what() const noexcept override { return m.c_str(); }
    };
}

namespace housepi::msg
{
    /// Show a message to stdout and flush.
    void info(const std::string &txt) noexcept;
    /// Display a warning message if not already displayed before.
    void warn(const std::string &txt);
    /// Wheter text represent a previously unseen warning.
    bool is_new_warning(const std::string &txt) noexcept;
    /// Mark given text as an already seen warning.
    void add_new_warning(const std::string &txt) noexcept;
    /// Throw an error message.
    [[noreturn]] void stop(const std::string &txt);

    /// cli inline-markup builders. Wrap identifier-safe tokens only.
    inline std::string arg(const std::string &x) { return "{.arg " + x + "}"; }
    inline std::string fn(const std::string &x) { return "{.fn " + x + "}"; }
    inline std::string val(const std::string &x) { return "{.val " + x + "}"; }
    inline std::string field(const std::string &x) { return "{.field " + x + "}"; }
    inline std::string code(const std::string &x) { return "{.code " + x + "}"; }
}
