/* housepi - Epidemiological Analysis of Household Outbreaks.
 *
 * Copyright © 2024-2026 Institut Pasteur (Paris), CNRS and INSERM.
 * Distributed under the terms of the GNU General Public License (GPLv3).
 * See the LICENSE file for details.
 */

#pragma once

#include <algorithm>
#include <unordered_set>
#include <vector>

// UPDATE: these operators should be cleaned with a proper namespace...

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
    std::vector<T> vec = x;
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
    for (const T &item : v)
        if (!exclusions.count(item))
            result.push_back(item);
    return result;
}

/// Whether the given element is contained in the given vector.
template <typename T, typename U>
bool contains(const T &vec, U elem)
{
    return std::find(vec.begin(), vec.end(), elem) != vec.end();
}

/// Appends elements from another container that are not already in the first.
template <typename T>
void union_append(T &vec, const T &other) { vec += (other - vec); }
