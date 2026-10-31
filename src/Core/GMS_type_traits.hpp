
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

#ifndef __GMS_TYPE_TRAITS_HPP__
#define __GMS_TYPE_TRAITS_HPP__ 311020261246

#include <boost/multi_array.hpp>
#include <vector>
#include <valarray>
#include <type_traits>

namespace gms
{

namespace common
{

template <typename T, class Container_t>
struct Is_Container_std_vector : public std::integral_constant<bool, false> {};

template <typename T>
struct Is_Container_std_vector<T, std::vector<T>>
    : public std::integral_constant<bool, true> {};

template <typename T, class Container_t>
struct Is_Container_std_valarray 
    : public std::integral_constant<bool, false> {};

template <typename T>
struct Is_Container_std_valarray<T, std::valarray<T>>
    : public std::integral_constant<bool, true> {};

template <typename T, class Container_t>
struct Is_Container_boost_multi_array1d
    : public std::integral_constant<bool, false> {};

template <typename T>
struct Is_Container_boost_multi_array1d<T, boost::multi_array<T, 1>>
    : public std::integral_constant<bool, true> {};

}

}
#endif /*__GMS_TYPE_TRAITS_HPP__ */