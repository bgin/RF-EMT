
/* Copyright (C) Bernard Gingold, 2020-2026 
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef __GMS_TYPE_ALIASES_HPP__
#define __GMS_TYPE_ALIASES_HPP__ 021120260540

#include <type_traits>

namespace gms
{

namespace common
{

template <bool B, typename T = void>
using ENABLE_IF = typename std::enable_if<B, T>::type;

template <typename T> using IS_CLASS = std::is_class<T>;

template <typename T> using IS_FLOATING_POINT = std::is_floating_point<T>;

template <typename T> using IS_INTEGRAL = std::is_integral<T>;

template <typename T> using IS_OBJECT = std::is_object<T>;

template <typename T> using IS_UNION = std::is_union<T>;

template <bool B, class T, class F>
using CONDITIONAL = std::conditional<B, T, F>;


template <typename T>
auto Wrf_Is_floating_point() noexcept(true)
    -> ENABLE_IF<IS_FLOATING_POINT<T>::value, constexpr bool> {

  return (IS_FLOATING_POINT<T>::value);
}

template <typename T>
auto Is_class() noexcept(true)
    -> ENABLE_IF<IS_CLASS<T>::value, constexpr bool> {

  return (IS_CLASS<T>::value);
}

template <typename T>
auto Is_object() noexcept(true)
    -> ENABLE_IF<IS_OBJECT<T>::value, constexpr bool> {

  return (IS_OBJECT<T>::value);
}

template <typename T>
auto Is_union() noexcept(true)
    -> ENABLE_IF<IS_UNION<T>::value, constexpr bool> {

  return (IS_UNION<T>::value);
}

}

}

#endif /*__GMS_TYPE_ALIASES_HPP__*/