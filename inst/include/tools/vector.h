/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <algorithm>
#include <unordered_set>
#include <vector>

namespace housepi::vect
{
    /// Append `y` to `x`.
    template <typename T>
    void append(std::vector<T> &x, const std::vector<T> &y)
    {
        x.reserve(x.size() + y.size());
        x.insert(x.end(), y.begin(), y.end());
    }

    /// Return `y` appended to `x`.
    template <typename T>
    std::vector<T> concatenate(const std::vector<T> &x, const std::vector<T> &y)
    {
        std::vector<T> out = x;
        append(out, y);
        return out;
    }

    /// Remove from `v` every element of `w`.
    template <typename T>
    void remove(std::vector<T> &v, const std::vector<T> &w)
    {
        const std::unordered_set<T> exclusions(w.begin(), w.end());
        v.erase(std::remove_if(v.begin(), v.end(),
                               [&exclusions](const T &item) { return exclusions.count(item); }),
                v.end());
    }

    /// Return `v` without the elements of `w`.
    template <typename T>
    std::vector<T> difference(const std::vector<T> &v, const std::vector<T> &w)
    {
        std::vector<T> out = v;
        remove(out, w);
        return out;
    }

    /// Whether `vec` contains `elem`.
    template <typename T, typename U>
    bool contains(const T &vec, U elem)
    {
        return std::find(vec.begin(), vec.end(), elem) != vec.end();
    }

    /// Append to `v` the elements of `other` absent from `v`.
    template <typename T>
    void append_unique(std::vector<T> &v, const std::vector<T> &other)
    {
        append(v, difference(other, v));
    }
}
