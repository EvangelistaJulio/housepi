/* housepi - Epidemiological Analysis of Household Outbreaks.
 * Copyright 2024-2026 Institut Pasteur (Paris), France.
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <algorithm>
#include <vector>

namespace housepi
{
    using matrix = std::vector<std::vector<double>>; /// Dense table of doubles, as equal-length vectors.
}

namespace housepi::vect
{
    /// Append `y` to `x`.
    template <typename T>
    void append(std::vector<T> &x, const std::vector<T> &y)
    {
        x.insert(x.end(), y.begin(), y.end());
    }

    /// Return `y` appended to `x`.
    template <typename T>
    std::vector<T> concatenate(std::vector<T> x, const std::vector<T> &y)
    {
        append(x, y);
        return x;
    }

    /// Whether `vec` contains `elem`.
    template <typename T, typename U>
    bool contains(const T &vec, const U &elem)
    {
        return std::find(vec.begin(), vec.end(), elem) != vec.end();
    }

    /// Sort `v` in increasing order, dropping duplicates.
    template <typename T>
    void sort_unique(std::vector<T> &v)
    {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }

    /// Remove from `v` every element of `w`.
    template <typename T>
    void remove(std::vector<T> &v, std::vector<T> w)
    {
        sort_unique(w);
        v.erase(std::remove_if(v.begin(), v.end(),
                               [&w](const T &x)
                               { return std::binary_search(w.begin(), w.end(), x); }),
                v.end());
    }

    /// Return `v` without the elements of `w`.
    template <typename T>
    std::vector<T> difference(std::vector<T> v, const std::vector<T> &w)
    {
        remove(v, w);
        return v;
    }

    /// Append to `v` the elements of `other` absent from `v`.
    template <typename T>
    void append_unique(std::vector<T> &v, const std::vector<T> &other)
    {
        append(v, difference(other, v));
    }
}
