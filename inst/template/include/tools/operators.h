/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#ifndef HOUSEPI__TOOLS__OPERATORS__H
#define HOUSEPI__TOOLS__OPERATORS__H

#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include "tools/types.h"

/// (% overload) Check if container containes given value.
template <typename Container, typename T>
bool operator%(const Container &c, const T &x)
{
    return std::find(c.begin(), c.end(), x) != c.end();
}

/// (+ overload) Concatenate containers.
template <typename T>
std::vector<T> operator+(const std::vector<T> &x, const std::vector<T> &y)
{
    std::vector<T> vec;
    vec.reserve(x.size() + y.size());
    vec.insert(vec.end(), x.begin(), x.end());
    vec.insert(vec.end(), y.begin(), y.end());
    return vec;
}

/// (+ overload) Concatenate containers.
template <typename T>
std::vector<T> &operator+=(std::vector<T> &x, const std::vector<T> &y)
{
    x.reserve(x.size() + y.size());
    x.insert(x.end(), y.begin(), y.end());
    return x;
}

/// (- overload) Remove all elements of first containers contained in second.
template <typename T>
std::vector<T> operator-(std::vector<T> v, const std::vector<T> &w)
{
    auto new_end = std::remove_if(v.begin(), v.end(),
                                  [&w](const T &i)
                                  { return w % i; });
    v.erase(new_end, v.end());
    return std::move(v);
}

namespace hsp
{
    /// Safely access an element in a u_dict by key.
    /// If the key does not exist, an exception is thrown with the key name.
    /// UPDATE: implement heterogeneous lookup to get rid of temporary string (C++20)
    template <typename T>
    inline T &safe_at(dict<T> &m, const std::string &k, const char *what = "key")
    {
        const auto it = m.find(k);
        if (it != m.end())
            return it->second;
        throw std::out_of_range(std::string("Unknown ") + what + " '" + k + "'.");
    }
    template <typename T>
    inline const T &safe_at(const dict<T> &m, const std::string &k, const char *what = "key")
    {
        const auto it = m.find(k);
        if (it != m.end())
            return it->second;
        throw std::out_of_range(std::string("Unknown ") + what + " '" + k + "'.");
    }

    /// Extract keys from dictionary.
    template <typename T>
    std::vector<std::string> extract_keys(const dict<T> &input)
    {
        std::vector<std::string> entries;
        entries.reserve(input.size());
        for (const auto &[key, val] : input)
            entries.push_back(key);
        return entries;
    }
}

#endif // HOUSEPI__TOOLS__OPERATORS__H
