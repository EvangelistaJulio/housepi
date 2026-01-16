/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__KEYS__H
#define HOUSEPI__TOOLS__KEYS__H

#include <stdexcept>
#include <string_view>
#include <vector>

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

        void init(std::vector<std::string_view> keys) { string_keys = std::move(keys); }

        const std::vector<std::string_view> &get_all_keys() const { return string_keys; }
        std::string_view get_key_name(T key) const { return get_key_name(static_cast<size_t>(key)); }
        std::string_view get_key_name(size_t index) const
        {
            if (index < string_keys.size())
                return string_keys[index];
            throw std::out_of_range("Key index out of range.");
        }

    private:
        std::vector<std::string_view> string_keys;
    };

    /// Parse comma-separated enum names from stringified VA_ARGS.
    std::vector<std::string_view> helper_parse_names(std::string_view input);
}

enum class AugmKeys : size_t;
enum class DataKeys : size_t;
enum class GroupRateKeys : size_t;
enum class ParamKeys : size_t;

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

#endif // HOUSEPI__TOOLS__KEYS__H