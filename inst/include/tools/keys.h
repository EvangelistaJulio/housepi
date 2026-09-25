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
    struct KeyHolder
    {
    public:
        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Assign the vector of model declared keys.
        void init(const std::vector<std::string> &keys);
        /// Drop all optional definitions.
        void clear_optional();
        /// Add a new definition.
        void add(std::string_view name);

        //////////////////////
        //////  GETTER  //////
        //////////////////////

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

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        /// Whether no definition has been provided.
        bool empty() const noexcept { return size() == 0; };
        /// Wheter a given definition is registred.
        bool count(std::string_view name) const;

    private:
        size_t size_ = 0;                        /// Store registry size
        std::vector<std::string> declared_ = {}; /// Model declared definitions
        std::vector<std::string> optional_ = {}; /// User optional definitions
    };

    struct KeyModule
    {
    public:
        /// Initialize key-holder with model declared keys.
        void init_key_holder(const std::vector<std::string> &keys) { khold_.init(keys); }

        /// Reference to key-holder.
        const KeyHolder &khold() const { return khold_; }
        /// Call key-holder index-to-string converter.
        std::string khold(size_t index) const { return khold_.as_string(index); }
        std::vector<std::string> khold(std::vector<size_t> indexes) const { return khold_.as_string(indexes); }
        /// Call key-holder string-to-index converter.
        size_t khold(std::string_view name) const { return khold_.from_string(name); }
        std::vector<size_t> khold(const std::vector<std::string> &names) const { return khold_.from_string(names); }

    protected:
        KeyHolder khold_; /// definitions key-holder
    };

    struct KeyModuleWithRequests : public KeyModule
    {
        /// Request discrete definition.
        void request_discrete(std::vector<size_t> vec);
        /// Request positive definition.
        void request_positive(std::vector<size_t> vec);

        /// Whether a key was requested discrete.
        bool is_required_discrete(size_t p) const noexcept;
        /// Whether a key was requested positive.
        bool is_required_positive(size_t p) const noexcept;
        /// Names of the keys requested discrete.
        std::vector<std::string> discrete_keys() const { return khold(discr_); }
        /// Names of the keys requested positive.
        std::vector<std::string> positive_keys() const { return khold(posit_); }

        /// Copy the requests registered on a module over the same keys.
        void copy_requests(const KeyModuleWithRequests &other)
        {
            discr_ = other.discr_;
            posit_ = other.posit_;
        }

        /// Throw if a definition does not satisfy the requests registred on its key.
        template <class Def>
        void assert_requests(size_t p, const Def &def) const;

    private:
        std::vector<size_t> discr_ = {}; /// Keys requested discrete
        std::vector<size_t> posit_ = {}; /// Keys requested positive
    };
}
