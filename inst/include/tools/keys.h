/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace hsp::key
{
    struct Registry
    {
    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        Registry(const std::vector<std::string> &keys) : declared_(keys) {}

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        /// Drop all optional definitions.
        void clear() { optional_.clear(); }
        /// Add a new definition.
        void add(std::string_view name)
        {
            if (!count(name))
                optional_.emplace_back(name);
        }

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        /// Total number of definition names.
        size_t size() const noexcept { return declared_.size() + optional_.size(); }

        /// Return all definition names.
        std::vector<std::string> get_all() const;

        /// Convert index to definition name.
        std::string as_string(size_t index) const;
        std::vector<std::string> as_string(std::vector<size_t> indexes);

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
        std::vector<std::string> declared_ = {}; /// Model declared definitions
        std::vector<std::string> optional_ = {}; /// User optional definitions
    };
}
