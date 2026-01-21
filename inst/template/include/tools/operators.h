/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <algorithm>
#include <unordered_set>
#include <vector>

/// (+= overload) Concatenate containers.
template <typename T>
std::vector<T> &operator+=(std::vector<T> &x, const std::vector<T> &y)
{
    x.reserve(x.size() + y.size());
    x.insert(x.end(), y.begin(), y.end());
    return x;
}

/// (+ overload) Concatenate containers.
template <typename T>
std::vector<T> operator+(const std::vector<T> &x, const std::vector<T> &y)
{
    auto vec = x;
    vec += y;
    return vec;
}

/// (- overload) Remove all elements of first containers contained in second.
template <typename T>
std::vector<T> operator-(const std::vector<T> &v, const std::vector<T> &w)
{
    const std::unordered_set<T> exclusions(w.begin(), w.end());
    std::vector<T> result;
    result.reserve(v.size());
    for (const auto &item : v)
        if (!exclusions.count(item))
            result.push_back(item);
    return result;
}
