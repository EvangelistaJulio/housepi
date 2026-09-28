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
    /// Finite set of entity type names displayed in errors.
    enum class ErrType : std::uint8_t
    {
        Empty,
        Numeric,
        Pointer,
        Function,
        Parameter,
        Latent,
        Observation,
        HouseholdObservation,
        GroupRisk,
        Data,
        Distribution,
        InfecProfile,
        Proposal
    };

    namespace str
    {
        /// Display word for error type `type`, e.g. "household observation".
        std::string format(ErrType type) noexcept;
    }
}

namespace housepi::msg
{
    /// Print message `txt` to stdout and flush.
    void info(const std::string &txt) noexcept;
    /// Raise warning `txt`.
    void warn(const std::string &txt);
    /// Throw error `txt`.
    [[noreturn]] void stop(const std::string &txt);

    /// Build cli inline markup.
    /// @note Wrap identifier-safe tokens only.
    inline std::string arg(const std::string &x) { return "{.arg " + x + "}"; }
    inline std::string fn(const std::string &x) { return "{.fn " + x + "}"; }
    inline std::string val(const std::string &x) { return "{.val " + x + "}"; }
    inline std::string field(const std::string &x) { return "{.field " + x + "}"; }
    inline std::string code(const std::string &x) { return "{.code " + x + "}"; }
    /// Tag a message line as an alert bullet (`!`) instead of the default hint.
    inline std::string alert(const std::string &x) { return "\x1F!" + x; }
    /// Set the cli pluralization count of a following `{?singular/plural}`.
    template <class T>
    inline std::string qty(const T &v) { return "{cli::qty(" + std::to_string(v.size()) + ")}"; }
    
    /// Append context `ctx` to the last sentence of message `what`.
    inline std::string join(std::string what, const std::string &ctx)
    {
        if (!what.empty() && what.back() == '.')
        {
            what.back() = ' ';
            return what + ctx + ".";
        }
        return what + "\n" + ctx + ".";
    }
}
