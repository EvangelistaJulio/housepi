/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace housepi
{
    /// Register definition names: model-declared ones first, then optional user ones.
    struct KeyHolder
    {
    public:
        //////  Setter  //////

        /// Assign model-declared keys `keys`, dropping optional ones.
        void init(const std::vector<std::string> &keys);
        /// Drop all optional definitions.
        void clear_optional();
        /// Add optional definition `name`, if absent.
        void add(std::string_view name);

        //////  Getter  //////

        /// Total number of definition names.
        size_t size() const noexcept { return size_; }

        /// Return all definition names.
        std::vector<std::string> all() const;

        /// Convert index to definition name.
        std::string as_string(size_t index) const;
        std::vector<std::string> as_string(std::vector<size_t> indexes) const;

        /// Convert definition name to index.
        size_t from_string(std::string_view name) const;
        std::vector<size_t> from_string(const std::vector<std::string> &names) const;

        //////  Flags  //////

        /// Whether no definition has been provided.
        bool empty() const noexcept { return size() == 0; };
        /// Whether definition `name` is registered.
        bool count(std::string_view name) const;

    private:
        size_t size_ = 0;                        /// Number of definitions.
        std::vector<std::string> declared_ = {}; /// Model-declared definitions.
        std::vector<std::string> optional_ = {}; /// Optional user definitions.
    };

    /// Module addressing its definitions through a key holder.
    struct KeyModule
    {
    public:
        /// Initialize the key holder with model-declared keys.
        void init_key_holder(const std::vector<std::string> &keys) { khold_.init(keys); }

        /// Key holder.
        const KeyHolder &khold() const { return khold_; }
        /// Convert indexes to definition names.
        std::string khold(size_t index) const { return khold_.as_string(index); }
        std::vector<std::string> khold(std::vector<size_t> indexes) const { return khold_.as_string(indexes); }
        /// Convert definition names to indexes.
        size_t khold(std::string_view name) const { return khold_.from_string(name); }
        std::vector<size_t> khold(const std::vector<std::string> &names) const { return khold_.from_string(names); }

    protected:
        KeyHolder khold_; /// Key holder of the definitions.
    };

    /// Key module recording discrete and positive requests per key.
    struct KeyModuleWithRequests : public KeyModule
    {
        /// Request discrete definitions `vec`.
        void request_discrete(std::vector<size_t> vec);
        /// Request positive definitions `vec`.
        void request_positive(std::vector<size_t> vec);

        /// Whether a key was requested discrete.
        bool is_required_discrete(size_t p) const noexcept;
        /// Whether a key was requested positive.
        bool is_required_positive(size_t p) const noexcept;
        /// Names of the keys requested discrete.
        std::vector<std::string> discrete_keys() const { return khold(discr_); }
        /// Names of the keys requested positive.
        std::vector<std::string> positive_keys() const { return khold(posit_); }

        /// Copy the requests registered on `other`, a module over the same keys.
        void copy_requests(const KeyModuleWithRequests &other)
        {
            discr_ = other.discr_;
            posit_ = other.posit_;
        }

        /// Throw if definition `def` does not satisfy the requests registered on key `p`.
        template <class Def>
        void assert_requests(size_t p, const Def &def) const;

    private:
        std::vector<size_t> discr_ = {}; /// Keys requested discrete.
        std::vector<size_t> posit_ = {}; /// Keys requested positive.
    };
}
