/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <cstdint>
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

    /// Finite set of entity type names displayed in errors.
    enum class ErrType : std::uint8_t
    {
        Empty,
        Numeric,
        Pointer,
        Function,
        Parameter,
        Augmented,
        Observation,
        HouseholdObservation,
        GroupRisk,
        Data,
        Distribution,
        InfecProfile,
        Proposal
    };

    /// Display word for an error type (e.g. "household observation").
    std::string err_type_name(ErrType type) noexcept;
}

namespace housepi::msg
{
    /// Show a message to stdout and flush.
    void info(const std::string &txt) noexcept;
    /// Display a warning message.
    void warn(const std::string &txt);
    /// Throw an error message.
    [[noreturn]] void stop(const std::string &txt);

    /// cli inline-markup builders. Wrap identifier-safe tokens only.
    inline std::string arg(const std::string &x) { return "{.arg " + x + "}"; }
    inline std::string fn(const std::string &x) { return "{.fn " + x + "}"; }
    inline std::string val(const std::string &x) { return "{.val " + x + "}"; }
    inline std::string field(const std::string &x) { return "{.field " + x + "}"; }
    inline std::string code(const std::string &x) { return "{.code " + x + "}"; }
    /// cli pluralization quantity: sets the count for a following `{?singular/plural}`.
    template <class T>
    inline std::string qty(const T &v) { return "{cli::qty(" + std::to_string(v.size()) + ")}"; }
}
