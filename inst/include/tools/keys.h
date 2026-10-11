/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "tools/nan.h"

namespace housepi
{
    /// Register definition names: model-declared ones first, then optional user ones.
    ///
    /// Modules derive from it:
    /// - read: public;
    /// - write: by the deriving module only.
    struct KeyHolder
    {
    public:
        //////  Constructor  //////

        /// Register model-declared keys `keys`.
        explicit KeyHolder(std::vector<std::string> keys)
            : keys_(std::move(keys)), n_declared_(keys_.size()) {}

        //////  Getter  //////

        /// Total number of definition names.
        size_t size() const noexcept { return keys_.size(); }

        /// Return all definition names.
        const std::vector<std::string> &all() const noexcept { return keys_; }

        /// Convert index to definition name.
        std::string consult(size_t index) const;
        std::vector<std::string> consult(const std::vector<size_t> &indexes) const;
        /// Convert definition name to index.
        /// - name not registered: error.
        size_t consult(std::string_view name) const;
        std::vector<size_t> consult(const std::vector<std::string> &names) const;
        /// Convert definition name to index.
        /// - name not registered: NA.
        size_t find(std::string_view name) const noexcept;
        std::vector<size_t> find(const std::vector<std::string> &names) const;

        //////  Flags  //////

        /// Whether no definition has been provided.
        bool empty() const noexcept { return size() == 0; };
        /// Whether definition `name` is registered.
        bool count(std::string_view name) const noexcept { return !na::is(find(name)); }

    protected:
        //////  Setter  //////

        /// Drop all optional definitions.
        void clear_optional();
        /// Add optional definition `name`, if absent.
        void add(std::string_view name);

    private:
        std::vector<std::string> keys_ = {}; /// Model-declared definitions, then optional ones.
        size_t n_declared_ = 0;              /// Number of model-declared definitions.
    };

    /// Key holder recording discrete, positive and boolean requests per key.
    struct KeyHolderWithRequests : public KeyHolder
    {
        using KeyHolder::KeyHolder;

        /// Request discrete definitions `vec`.
        void request_discrete(std::vector<size_t> vec);
        /// Request positive definitions `vec`.
        void request_positive(std::vector<size_t> vec);
        /// Request boolean definitions `vec`, valued 0 or 1.
        void request_boolean(std::vector<size_t> vec);

        /// Whether a key was requested discrete.
        bool is_required_discrete(size_t p) const noexcept;
        /// Whether a key was requested positive.
        bool is_required_positive(size_t p) const noexcept;
        /// Whether a key was requested boolean.
        bool is_required_boolean(size_t p) const noexcept;
        /// Names of the keys requested discrete.
        std::vector<std::string> discrete_keys() const { return consult(discr_); }

        /// Copy the requests registered on `other`, a key holder over the same keys.
        void copy_requests(const KeyHolderWithRequests &other)
        {
            discr_ = other.discr_;
            posit_ = other.posit_;
            boolean_ = other.boolean_;
        }

        /// Throw if definition `def` does not satisfy the requests registered on key `p`.
        template <class Def>
        void assert_requests(size_t p, const Def &def) const;

    private:
        std::vector<size_t> discr_ = {};   /// Keys requested discrete.
        std::vector<size_t> posit_ = {};   /// Keys requested positive.
        std::vector<size_t> boolean_ = {}; /// Keys requested boolean.
    };
}
