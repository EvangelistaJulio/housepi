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

namespace housepi::vect
{
    /// Appends the second container to the first.
    template <typename T>
    void append(std::vector<T> &x, const std::vector<T> &y)
    {
        x.reserve(x.size() + y.size());
        x.insert(x.end(), y.begin(), y.end());
    }

    /// Concatenation of two containers.
    template <typename T>
    std::vector<T> concatenate(const std::vector<T> &x, const std::vector<T> &y)
    {
        std::vector<T> out = x;
        append(out, y);
        return out;
    }

    /// Removes from the first container every element of the second.
    template <typename T>
    void remove(std::vector<T> &v, const std::vector<T> &w)
    {
        const std::unordered_set<T> exclusions(w.begin(), w.end());
        v.erase(std::remove_if(v.begin(), v.end(),
                               [&exclusions](const T &item) { return exclusions.count(item); }),
                v.end());
    }

    /// First container deprived of all the elements of the second.
    template <typename T>
    std::vector<T> difference(const std::vector<T> &v, const std::vector<T> &w)
    {
        std::vector<T> out = v;
        remove(out, w);
        return out;
    }

    /// Whether the given element is contained in the given vector.
    template <typename T, typename U>
    bool contains(const T &vec, U elem)
    {
        return std::find(vec.begin(), vec.end(), elem) != vec.end();
    }

    /// Appends elements from another container that are not already in the first.
    template <typename T>
    void append_unique(std::vector<T> &v, const std::vector<T> &other)
    {
        append(v, difference(other, v));
    }
}
