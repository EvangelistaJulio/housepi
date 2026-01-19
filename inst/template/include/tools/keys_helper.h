/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__KEYS_HELPER__H
#define HOUSEPI__TOOLS__KEYS_HELPER__H

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
        void add_key(std::string name) { optional_.push_back(name); }

        size_t size() const { return declared_.size() + optional_.size(); }
        std::vector<std::string> get_all_keys() const
        {
            std::vector<std::string> vec;
            vec.reserve(declared_.size() + optional_.size());
            vec.insert(vec.end(), declared_.begin(), declared_.end());
            vec.insert(vec.end(), optional_.begin(), optional_.end());
            return vec;
        }
        std::string as_string(T key) const { return as_string(static_cast<size_t>(key)); }
        std::string as_string(size_t index) const
        {
            if (index < declared_.size())
                return declared_[index];
            index -= declared_.size();
            if (index < optional_.size())
                return optional_[index];
            throw std::out_of_range("Key index out of range.");
        }
        // Convert from string to enum. Throw exception if not found.
        size_t from_string(std::string name) const
        {
            auto it = std::find(declared_.begin(), declared_.end(), name);
            if (it != declared_.end())
                return static_cast<size_t>(std::distance(declared_.begin(), it));
            auto it2 = std::find(optional_.begin(), optional_.end(), name);
            if (it2 != optional_.end())
                return declared_.size() + static_cast<size_t>(std::distance(optional_.begin(), it2));
            throw std::out_of_range("Key name '" + std::string(name) + "' not found.");
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
    std::vector<std::string> get_all_keys() { return Registry<T>::instance().get_all_keys(); }

    template <typename T>
    std::string as_string(T key) { return Registry<T>::instance().as_string(key); }

    template <typename T>
    std::string as_string(size_t index) { return Registry<T>::instance().as_string(index); }

    template <typename T>
    size_t from_string(std::string name) { return Registry<T>::instance().from_string(name); }

    /// Parse comma-separated enum names from stringified VA_ARGS.
    std::vector<std::string> helper_parse_names(std::string input);
}

/// Define enum keys and auto-generate registration.
#define REGISTER_HOUSEPI_KEYS(T, ...)                                    \
    enum class T : size_t                                                \
    {                                                                    \
        __VA_ARGS__                                                      \
    };                                                                   \
                                                                         \
    namespace                                                            \
    {                                                                    \
        struct T##_AutoRegister                                          \
        {                                                                \
            T##_AutoRegister()                                           \
            {                                                            \
                auto keys = hsp::key::helper_parse_names(#__VA_ARGS__);  \
                hsp::key::Registry<T>::instance().init(std::move(keys)); \
            }                                                            \
        };                                                               \
        static T##_AutoRegister T##_housepi_auto_register;               \
    }

#endif // HOUSEPI__TOOLS__KEYS_HELPER__H