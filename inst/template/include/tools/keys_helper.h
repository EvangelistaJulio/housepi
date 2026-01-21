/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

// UPDATE move back to string_view once Variable constructor is fixed
namespace hsp::key
{
    template <typename T>
    class Registry
    {
    public:
        static Registry &instance()
        {
            static Registry reg;
            return reg;
        }

        void clear() { optional_.clear(); }
        void init(std::vector<std::string> keys)
        {
            clear();
            declared_ = keys;
        }

        void add_key(std::string_view name) { optional_.emplace_back(name); }

        size_t size() const { return declared_.size() + optional_.size(); }

        std::vector<std::string_view> get_all_keys() const
        {
            std::vector<std::string_view> vec;
            vec.reserve(declared_.size() + optional_.size());
            for (auto &s : declared_)
                vec.push_back(s);
            for (auto &s : optional_)
                vec.push_back(s);
            return vec;
        }

        bool count(std::string_view name) const
        {
            return std::find(declared_.begin(), declared_.end(), name) != declared_.end() ||
                   std::find(optional_.begin(), optional_.end(), name) != optional_.end();
        }

        std::string as_string(T key) const { return as_string(static_cast<size_t>(key)); }
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
            for (const auto &s : indexes)
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
            for (const auto &s : names)
                output.push_back(from_string(s));
            return output;
        }

    private:
        std::vector<std::string> declared_;
        std::vector<std::string> optional_;
    };

    template <typename T>
    void clear() { Registry<T>::instance().clear(); }

    template <typename T>
    size_t size() { return Registry<T>::instance().size(); }

    template <typename T>
    std::vector<std::string_view> get_all_keys() { return Registry<T>::instance().get_all_keys(); }

    template <typename T>
    bool count(std::string_view name) { return Registry<T>::instance().count(name); }

    template <typename T>
    std::string as_string(T key) { return Registry<T>::instance().as_string(key); }

    template <typename T>
    std::string as_string(size_t index) { return Registry<T>::instance().as_string(index); }

    template <typename T>
    std::vector<std::string> as_string(std::vector<size_t> indexes) { return Registry<T>::instance().as_string(indexes); }

    template <typename T>
    size_t from_string(std::string_view name) { return Registry<T>::instance().from_string(name); }

    template <typename T>
    std::vector<size_t> from_string(const std::vector<std::string> &names) { return Registry<T>::instance().from_string(names); }

    /// Parse comma-separated enum names from stringified VA_ARGS.
    std::vector<std::string> helper_parse_names(std::string input);
}
