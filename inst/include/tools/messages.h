/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace housepi
{
    /// Finite set of entity type names displayed in errors.
    enum class ErrType : std::uint8_t
    {
        Empty,
        Numeric,
        Pointer,
        Function,
        Operator,
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
    /// Print message `txt`, with no line ending added.
    void info(const std::string &txt);
    /// Raise warning `txt`.
    void warn(const std::string &txt);
    /// Throw error `txt`.
    [[noreturn]] void stop(const std::string &txt);

    /// Store the warnings raised by the running thread.
    ///
    /// * each distinct warning -> raised once, by flush();
    /// * nested scopes -> the innermost holds;
    /// * scope left without flush(), on an error -> warnings dropped.
    struct WarningBuffer
    {
        WarningBuffer() noexcept;
        ~WarningBuffer();
        WarningBuffer(const WarningBuffer &) = delete;
        WarningBuffer &operator=(const WarningBuffer &) = delete;

        /// Hold warning `txt` in the innermost buffer of the calling thread.
        /// @return Whether a buffer holds it.
        static bool hold(const std::string &txt);
        /// Leave the scope, then raise each distinct warning once.
        void flush();

    private:
        std::vector<std::string> warnings_ = {}; /// Distinct warnings, in order of first occurrence.
        WarningBuffer *outer_;                   /// Buffer in scope before this one.
        bool is_open_ = true;                    /// Whether it holds the warnings of its thread.
    };

    /// Build cli inline markup.
    /// @note Wrap identifier-safe tokens only.
    inline std::string arg(const std::string &x) { return "{.arg " + x + "}"; }
    inline std::string fn(const std::string &x) { return "{.fn " + x + "}"; }
    inline std::string val(const std::string &x) { return "{.val " + x + "}"; }
    inline std::string field(const std::string &x) { return "{.field " + x + "}"; }
    inline std::string code(const std::string &x) { return "{.code " + x + "}"; }
    /// Show housepi call `call`, run on click.
    inline std::string run(const std::string &call) { return "{.run [" + call + "](housepi::" + call + ")}"; }
    /// Show `?topic`, opening the housepi help page on click.
    inline std::string help(const std::string &topic)
    {
        return "{.run [?" + topic + "](utils::help(\"" + topic + "\", \"housepi\"))}";
    }
    /// Tag a message line as an alert bullet (`!`) instead of the default hint.
    inline std::string alert(const std::string &x) { return "\x1F!" + x; }
    /// Set the cli pluralization count of a following `{?singular/plural}`.
    template <class T>
    inline std::string qty(const T &v) { return "{cli::qty(" + std::to_string(v.size()) + ")}"; }

    /// Append context `ctx` to the first line of message `what`, ahead of its hint lines.
    inline std::string join(std::string what, const std::string &ctx)
    {
        size_t end = what.find('\n');
        if (end == std::string::npos)
            end = what.size();
        if (end > 0 && what[end - 1] == '.')
            what.erase(--end, 1);
        return what.insert(end, " " + ctx + ".");
    }
}
