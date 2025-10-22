/*
 *        __
 *       / /_  ____  __  __ ____ __  ____  __
 *      / __ \/ __ \/ / / / ___/ _ \/ __ \/ /
 *     / / / / /_/ / /_/ (__  )  __/ /_/ / /
 *    /_/ /_/\____/\____/____/\___/ ____/_/
 *                               /_/
 *
 * Copyright (C) 2025 Francesco Camaglia, Institut Pasteur
 *
 * tools/operators.h
 *
 */

#ifndef HOUSEPI__TOOLS__OPERATORS__H
#define HOUSEPI__TOOLS__OPERATORS__H

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
    return v;
}

/// Safely access an element in a u_dict by key.
/// If the key does not exist, an exception is thrown with the key name.
template <typename T>
inline T &safe_at(u_dict<T> &m, const std::string &k, const char *what = "key")
{
    const auto it = m.find(k);
    if (it != m.end())
        return it->second;
    throw std::out_of_range(std::string("Unknown ") + what + " '" + k + "'.");
}
template <typename T>
inline const T &safe_at(const u_dict<T> &m, const std::string &k, const char *what = "key")
{
    const auto it = m.find(k);
    if (it != m.end())
        return it->second;
    throw std::out_of_range(std::string("Unknown ") + what + " '" + k + "'.");
}

#endif // HOUSEPI__TOOLS__OPERATORS__H
