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
    /// UPDATE: enforce string_view once Variable constructor is fixed
    class Registry
    {
    public:
        ///////////////////////////
        //////  CONSTRUCTOR  //////
        ///////////////////////////

        Registry(const std::vector<std::string> &keys) : declared_(keys) { clear(); }

        //////////////////////
        //////  SETTER  //////
        //////////////////////

        void clear() { optional_.clear(); }

        void add(std::string_view name)
        {
            if (!count(name))
                optional_.emplace_back(name);
        }

        //////////////////////
        //////  GETTER  //////
        //////////////////////

        size_t size() const noexcept { return declared_.size() + optional_.size(); }

        std::vector<std::string> get_all() const;

        std::string as_string(size_t index) const;
        std::vector<std::string> as_string(std::vector<size_t> indexes);

        size_t from_string(std::string_view name) const;
        std::vector<size_t> from_string(const std::vector<std::string> &names) const;

        /////////////////////
        //////  FLAGS  //////
        /////////////////////

        bool count(std::string_view name) const;

    private:
        std::vector<std::string> declared_;
        std::vector<std::string> optional_;
    };
}
