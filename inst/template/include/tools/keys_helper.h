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

// FLAG move back to string_view once Variable constructor is fixed
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

        void init(std::vector<std::string> keys) { string_keys_ = keys; }

        size_t size() const { return string_keys_.size(); }
        const std::vector<std::string> &get_all_keys() const { return string_keys_; }
        std::string to_string(T key) const { return to_string(static_cast<size_t>(key)); }
        std::string to_string(size_t index) const
        {
            if (index < string_keys_.size())
                return string_keys_[index];
            throw std::out_of_range("Key index out of range.");
        }
        // Convert from string to enum. Throw exception if not found.
        T from_string(std::string name) const
        {
            auto it = std::find(string_keys_.begin(), string_keys_.end(), name);
            if (it != string_keys_.end())
                return static_cast<T>(std::distance(string_keys_.begin(), it));
            throw std::out_of_range("Key name '" + std::string(name) + "' not found.");
        }
        // Check if a key name exists.
        bool contains(std::string name) const
        {
            return std::find(string_keys_.begin(), string_keys_.end(), name) != string_keys_.end();
        }

    private:
        std::vector<std::string> string_keys_;
    };

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