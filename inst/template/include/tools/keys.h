/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace hsp::key
{
    /// UPDATE: enforce string_view once Variable constructor is fixed
    class Registry
    {
    public:
        Registry(const std::vector<std::string> &keys) : declared_(keys) { clear(); }

        void clear() { optional_.clear(); }

        void add(std::string_view name)
        {
            if (std::find(optional_.begin(), optional_.end(), name) == optional_.end())
                optional_.emplace_back(name);
        }

        size_t size() const noexcept { return declared_.size() + optional_.size(); }

        std::vector<std::string> get_all() const
        {
            std::vector<std::string> vec;
            vec.reserve(declared_.size() + optional_.size());
            for (const std::string &s : declared_)
                vec.push_back(s);
            for (const std::string &s : optional_)
                vec.push_back(s);
            return vec;
        }

        bool count(std::string_view name) const
        {
            return std::find(declared_.begin(), declared_.end(), name) != declared_.end() ||
                   std::find(optional_.begin(), optional_.end(), name) != optional_.end();
        }

        std::string as_string(size_t index) const
        {
            if (index < declared_.size())
                return declared_[index];
            index -= declared_.size();
            if (index < optional_.size())
                return optional_[index];
            throw std::out_of_range("Registration failed: key index out of range.");
        }
        std::vector<std::string> as_string(std::vector<size_t> indexes)
        {
            std::vector<std::string> output;
            output.reserve(indexes.size());
            for (size_t s : indexes)
                output.push_back(as_string(s));
            return output;
        }

        size_t from_string(std::string_view name) const
        {
            auto it = std::find(declared_.begin(), declared_.end(), name);
            if (it != declared_.end())
                return static_cast<size_t>(std::distance(declared_.begin(), it));
            auto it2 = std::find(optional_.begin(), optional_.end(), name);
            if (it2 != optional_.end())
                return declared_.size() + static_cast<size_t>(std::distance(optional_.begin(), it2));
            throw std::out_of_range("Registration failed: key name '" + std::string(name) + "' was not found.");
        }

        std::vector<size_t> from_string(const std::vector<std::string> &names) const
        {
            std::vector<size_t> output;
            output.reserve(names.size());
            for (const std::string &s : names)
                output.push_back(from_string(s));
            return output;
        }

    private:
        std::vector<std::string> declared_;
        std::vector<std::string> optional_ = {};
    };

    /// Parse comma-separated enum names from stringified variadic arguments.
    std::vector<std::string> helper_parse_names(std::string input);
}
